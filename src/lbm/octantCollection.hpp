    #ifndef OCTANT_COLLECTION_HPP
#define OCTANT_COLLECTION_HPP

#include "octantCollection.h"

#include <p8est_extended.h>

#include "../mpi/mpiManager.h"


namespace ALBRT {
    
template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
int OctantCollection<T, Lattice, nx, ny, nz, nghost>::numLocalOctants = 0;

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
OctantCollection<T, Lattice, nx, ny, nz, nghost>::
OctantCollection(int numOctantsX, int numOctantsY, int numOctantsZ, T octantPhysXLength, int maxLevel,
                bool periodicX, bool periodicY, bool periodicZ)
:   maxLevel(maxLevel)
{
    baseOctantPhysLength[0] = octantPhysXLength;
    baseOctantPhysLength[1] = octantPhysXLength * (S)ny / (S)nx;
    baseOctantPhysLength[2] = octantPhysXLength * (S)nz / (S)nx;

    connectivity = p8est_connectivity_new_brick(numOctantsX, numOctantsY, numOctantsZ, 
                                                (int) periodicX, (int) periodicY, (int) periodicZ);

    forest = p8est_new_ext(MPIManager::mpiCommWorld, connectivity, 0, 0, 0, sizeof(OctantInternalInfo<T>), NULL, (void *) this);

    p8est_partition(forest, 1, NULL);
    p8est_balance(forest, P8EST_CONNECT_CORNER, NULL);
    p8est_partition(forest, 1, NULL);

    p8est_iterate(forest, NULL, NULL, initOctantIndices, NULL, NULL, NULL); //need to initialize octant information prior to ghost exchange

    forestGhost = p8est_ghost_new(forest, P8EST_CONNECT_FULL);
    ghostOctantInfos.reserve(forestGhost->ghosts.elem_count);
    p8est_ghost_exchange_data(forest, forestGhost, (void*) ghostOctantInfos.data());

    p8est_iterate(forest, forestGhost, (void*) ghostOctantInfos.data(), NULL, initOctantFaceNeighbors, initOctantEdgeNeighbors, initOctantCornerNeighbors);

    for (int rankId = 0; rankId < MPIManager::numRanks; rankId++) {
        if (MPIManager::rank == rankId) {
            for (OctantInternalInfo<T>& info : localOctantInfos) {
                printf("Octant ID: %d, proc: %d, level %d\n", info.index, info.proc, info.level);
                for (int iF = 0; iF < FaceCommLattice::q; iF++) {
                    printf("Face neighbor %d: %s", iF, neighborTypeStrings[info.faceNeighbors[iF].type]);
                    if (info.faceNeighbors[iF].type == MATCH) {
                        printf(", neighbor proc %d, index: %d", info.faceNeighbors[iF].neighbor.match.proc, info.faceNeighbors[iF].neighbor.match.index);
                    }
                    printf("\n");
                }
            }
        }
        sc_MPI_Barrier(MPIManager::mpiCommWorld);
    }
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
OctantCollection<T, Lattice, nx, ny, nz, nghost>::~OctantCollection()
{
    p8est_destroy(forest);
    p8est_ghost_destroy(forestGhost);
    p8est_connectivity_destroy(connectivity);
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
void OctantCollection<T, Lattice, nx, ny, nz, nghost>::initOctantIndices(p8est_iter_volume_info_t *info, void *user_data)
{
    OctantCollection<T,Lattice,nx,ny,nz,nghost> * octCollection = 
        (OctantCollection<T,Lattice,nx,ny,nz,nghost> *) info->p4est->user_pointer;
    
    OctantInternalInfo<T> * octInfo = (OctantInternalInfo<T> *) info->quad->p.user_data;
    octInfo->proc = MPIManager::rank;
    octInfo->index = numLocalOctants;
    octInfo->level = info->quad->level;

    // Convert p4est internal qcoord system to vertex. Then scale using our physical octant parameters.
    double vxyz[3];
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, info->quad->x, info->quad->y, info->quad->z, vxyz);
    
    octInfo->origin[0] = octCollection->baseOctantPhysLength[0] * (T) vxyz[0];
    octInfo->origin[1] = octCollection->baseOctantPhysLength[1] * (T) vxyz[1];
    octInfo->origin[2] = octCollection->baseOctantPhysLength[2] * (T) vxyz[2];

    octInfo->physLength[0] = octCollection->baseOctantPhysLength[0] * std::pow(T(0.5), octInfo->level);
    octInfo->physLength[1] = octCollection->baseOctantPhysLength[1] * std::pow(T(0.5), octInfo->level);
    octInfo->physLength[2] = octCollection->baseOctantPhysLength[2] * std::pow(T(0.5), octInfo->level);
    
    octCollection->localOctantInfos.push_back(*octInfo);

    numLocalOctants++;
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
void OctantCollection<T, Lattice, nx, ny, nz, nghost>::initOctantFaceNeighbors(p8est_iter_face_info_t *info, void *user_data)
{
    OctantCollection<T,Lattice,nx,ny,nz,nghost> * octCollection = 
        (OctantCollection<T,Lattice,nx,ny,nz,nghost> *) info->p4est->user_pointer;

    OctantInternalInfo<T> * ghostInfos = (OctantInternalInfo<T> *) user_data;

    sc_array_t *sides = &(info->sides);
    p8est_iter_face_side_t *side[2];

    int numSides = sides->elem_count;

    for (int i = 0; i < numSides; i++) {
        side[i] = p8est_iter_fside_array_index_int(sides, i);
    }

    for (int i = 0; i < numSides; i++) {
        int face = side[i]->face;

        if (numSides == 1) { //process boundary face
            if (!side[i]->is.full.is_ghost) { //if this is not a ghost face (it shouldn't be)
                NeighborRelationship relationship;
                relationship.type = BOUNDARY;
                OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                info->faceNeighbors[face] = relationship;

                octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
            }
        }
        else if (side[i]->is_hanging) { //this side is split up into 4 quadrants
         
            for (int j = 0; j < numSides; j++) {
                if (i == j) continue;

                int compareFace = side[j]->face;
                if (face == FaceCommLattice::opp[compareFace]) {
                    assert(side[j]->is_hanging == false && "Only one side of a face interface should be hanging.");

                    NeighborRelationship relationship;
                    relationship.type = COARSE; //because side[i] is hanging, the type must be coarse.
                    
                    if (side[j]->is.full.is_ghost) {
                        OctantInternalInfo<T> * info = &ghostInfos[side[j]->is.full.quadid];
                        relationship.neighbor.coarse.info.proc = info->proc;
                        relationship.neighbor.coarse.info.index = info->index;
                    }
                    else {
                        OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[j]->is.full.quad->p.user_data;
                        relationship.neighbor.coarse.info.proc = info->proc;
                        relationship.neighbor.coarse.info.index = info->index;
                    }

                    for (int hangingQuad = 0; hangingQuad < 4; hangingQuad++) { //only operate on local quadrants
                        if (!side[i]->is.hanging.is_ghost[hangingQuad]) {
                            OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.hanging.quad[hangingQuad]->p.user_data;
                            relationship.neighbor.coarse.index = hangingQuad;
                            info->faceNeighbors[face] = relationship;

                            //update the vector as well.
                            octCollection->localOctantInfos.at(info->index) = *info;
                        }
                    }
                }
            }
        }
        else if (!side[i]->is.full.is_ghost) { //only work on local full faces
            for (int j = 0; j < numSides; j++) {
                if (i == j) continue;

                int compareFace = side[j]->face;
                if (face == FaceCommLattice::opp[compareFace]) {
                    if (side[j]->is_hanging) { //coarse octant neighbors 4 fine octants
                        NeighborRelationship relationship;
                        relationship.type = FINE;

                        for (int hangingQuad = 0; hangingQuad < 4; hangingQuad++) {
                            OctantInternalInfo<T> * info;
                            if (side[j]->is.hanging.is_ghost[hangingQuad])
                                info = (OctantInternalInfo<T> *) &ghostInfos[side[j]->is.hanging.quadid[hangingQuad]];
                            else
                                info = (OctantInternalInfo<T> *) side[j]->is.hanging.quad[hangingQuad]->p.user_data;

                            relationship.neighbor.fine[hangingQuad].proc = info->proc;
                            relationship.neighbor.fine[hangingQuad].index = info->index;

                        }

                        OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                        info->faceNeighbors[face] = relationship;

                        octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
                    }
                    else { //match relationship
                        NeighborRelationship relationship;
                        relationship.type = MATCH;

                        //get the matching neighbor octant
                        OctantInternalInfo<T> * info;
                        if (side[j]->is.full.is_ghost)
                            info = (OctantInternalInfo<T> *) &ghostInfos[side[j]->is.full.quadid];
                        else 
                            info = (OctantInternalInfo<T> *) side[j]->is.full.quad->p.user_data;

                        relationship.neighbor.match.proc = info->proc;
                        relationship.neighbor.match.index = info->index;

                        //reset info to the current octant
                        info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                        info->faceNeighbors[face] = relationship;

                        octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
                    }
                }
            }
        }
    }
        
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
void OctantCollection<T, Lattice, nx, ny, nz, nghost>::initOctantEdgeNeighbors(p8est_iter_edge_info_t *info, void *user_data)
{
    OctantCollection<T,Lattice,nx,ny,nz,nghost> * octCollection = 
        (OctantCollection<T,Lattice,nx,ny,nz,nghost> *) info->p4est->user_pointer;

    OctantInternalInfo<T> * ghostInfos = (OctantInternalInfo<T> *) user_data;

    sc_array_t *sides = &(info->sides);
    p8est_iter_edge_side_t *side[4];

    int numSides = sides->elem_count;

    for (int i = 0; i < numSides; i++) {
        side[i] = p8est_iter_eside_array_index_int(sides, i);
    }

    for (int i = 0; i < numSides; i++) {
        int edge = side[i]->edge;

        if (numSides == 1) { //process boundary face
            if (!side[i]->is.full.is_ghost) { //if this is not a ghost edge (it shouldn't be)
                NeighborRelationship relationship;
                relationship.type = BOUNDARY;
                OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                info->edgeNeighbors[edge] = relationship;

                octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
            }
        }
        else if (side[i]->is_hanging) {
            for (int j = 0; j < numSides; j++) {
                if (i == j) continue;

                int compareEdge = side[j]->edge;
                if (edge == EdgeCommLattice::opp[compareEdge]) {
                    assert(side[j]->is_hanging == false && "Only one side of an edge interface should be hanging.");

                    NeighborRelationship relationship;
                    relationship.type = COARSE; //because side[i] is hanging, the type must be coarse.
                    
                    if (side[j]->is.full.is_ghost) {
                        OctantInternalInfo<T> * info = &ghostInfos[side[j]->is.full.quadid];
                        relationship.neighbor.coarse.info.proc = info->proc;
                        relationship.neighbor.coarse.info.index = info->index;
                    }
                    else {
                        OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[j]->is.full.quad->p.user_data;
                        relationship.neighbor.coarse.info.proc = info->proc;
                        relationship.neighbor.coarse.info.index = info->index;
                    }

                    for (int hangingQuad = 0; hangingQuad < 2; hangingQuad++) { //only operate on local quadrants
                        if (!side[i]->is.hanging.is_ghost[hangingQuad]) {
                            OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.hanging.quad[hangingQuad]->p.user_data;
                            relationship.neighbor.coarse.index = hangingQuad;
                            info->edgeNeighbors[edge] = relationship;

                            //update the vector as well.
                            octCollection->localOctantInfos.at(info->index) = *info;
                        }
                    }
                }
            }
        }
        else if (!side[i]->is.full.is_ghost) { //only work on local full faces
            for (int j = 0; j < numSides; j++) {
                if (i == j) continue;

                int compareEdge = side[j]->edge;
                if (edge == EdgeCommLattice::opp[compareEdge]) {
                    if (side[j]->is_hanging) { //coarse octant neighbors 2 fine octants
                        NeighborRelationship relationship;
                        relationship.type = FINE;

                        for (int hangingQuad = 0; hangingQuad < 2; hangingQuad++) {
                            OctantInternalInfo<T> * info;
                            if (side[j]->is.hanging.is_ghost[hangingQuad])
                                info = (OctantInternalInfo<T> *) &ghostInfos[side[j]->is.hanging.quadid[hangingQuad]];
                            else
                                info = (OctantInternalInfo<T> *) side[j]->is.hanging.quad[hangingQuad]->p.user_data;

                            relationship.neighbor.fine[hangingQuad].proc = info->proc;
                            relationship.neighbor.fine[hangingQuad].index = info->index;

                        }

                        OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                        info->edgeNeighbors[edge] = relationship;

                        octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
                    }
                    else { //match relationship
                        NeighborRelationship relationship;
                        relationship.type = MATCH;

                        //get the matching neighbor octant
                        OctantInternalInfo<T> * info;
                        if (side[j]->is.full.is_ghost) 
                            info = (OctantInternalInfo<T> *) &ghostInfos[side[j]->is.full.quadid];
                        else 
                            info = (OctantInternalInfo<T> *) side[j]->is.full.quad->p.user_data;

                        relationship.neighbor.match.proc = info->proc;
                        relationship.neighbor.match.index = info->index;

                        //reset info to the current octant
                        info = (OctantInternalInfo<T> *) side[i]->is.full.quad->p.user_data;
                        info->edgeNeighbors[edge] = relationship;

                        octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
                    }
                }
            }
        }
    }        
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
void OctantCollection<T, Lattice, nx, ny, nz, nghost>::initOctantCornerNeighbors(p8est_iter_corner_info_t *info, void *user_data)
{
    OctantCollection<T,Lattice,nx,ny,nz,nghost> * octCollection = 
        (OctantCollection<T,Lattice,nx,ny,nz,nghost> *) info->p4est->user_pointer;

    OctantInternalInfo<T> * ghostInfos = (OctantInternalInfo<T> *) user_data;

    sc_array_t *sides = &(info->sides);
    p8est_iter_corner_side_t *side[8];

    int numSides = sides->elem_count;

    for (int i = 0; i < numSides; i++) {
        side[i] = p8est_iter_cside_array_index_int(sides, i);
    }

    for (int i = 0; i < numSides; i++) {
        int corner = side[i]->corner;

        if (numSides == 1) { //process boundary face
            if (!side[i]->is_ghost) { //if this is not a ghost corner (it shouldn't be)
                NeighborRelationship relationship;
                relationship.type = BOUNDARY;
                OctantInternalInfo<T> * info = (OctantInternalInfo<T> *) side[i]->quad->p.user_data;
                info->cornerNeighbors[corner] = relationship;

                octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
            }
        }
        else if (!side[i]->is_ghost) { //only work on local quads
            for (int j = 0; j < numSides; j++) {
                if (i == j) continue;

                int compareCorner = side[j]->corner;
                if (corner == CornerCommLattice::opp[compareCorner]) {
                    
                    NeighborRelationship relationship;

                    //get the matching neighbor octant
                    OctantInternalInfo<T> * info;
                    if (side[j]->is_ghost) 
                        info = (OctantInternalInfo<T> *) &ghostInfos[side[j]->quadid];
                    else 
                        info = (OctantInternalInfo<T> *) side[j]->quad->p.user_data;

                    int destLevel = info->level;                    

                    int proc = info->proc;
                    int index = info->index;

                    //reset info to the current octant
                    info = (OctantInternalInfo<T> *) side[i]->quad->p.user_data;
                    int level = info->level;

                    relationship.type = level == destLevel ? MATCH : (level < destLevel ? FINE : COARSE);

                    switch (relationship.type) {
                        case MATCH:
                            relationship.neighbor.match.index = index;
                            relationship.neighbor.match.proc = proc;
                            break;
                        case FINE:
                            relationship.neighbor.fine[0].index = index;
                            relationship.neighbor.fine[0].proc = proc;
                            break;
                        case COARSE:
                            relationship.neighbor.coarse.info.index = index;
                            relationship.neighbor.coarse.info.proc = proc;
                            relationship.neighbor.coarse.index = 0;
                            break;
                    }

                    info->cornerNeighbors[corner] = relationship;

                    octCollection->localOctantInfos.at(info->index) = *info; //update the collection vector as well
                }
            }
        }
    }        
}

} //end namespace

#endif