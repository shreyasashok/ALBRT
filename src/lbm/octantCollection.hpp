#ifndef OCTANT_COLLECTION_HPP
#define OCTANT_COLLECTION_HPP

#include "octantCollection.h"

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

    forest = p8est_new(MPIManager::mpiCommWorld, connectivity, sizeof(OctantInternalInfo<T>), NULL, NULL);
    forest->user_pointer = (void *) this;

    p8est_partition(forest, 1, NULL);
    p8est_balance(forest, P8EST_CONNECT_CORNER, NULL);
    p8est_partition(forest, 1, NULL);

    p8est_iterate(forest, NULL, NULL, initOctantIndices, NULL, NULL, NULL);

}

template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
OctantCollection<T, Lattice, nx, ny, nz, nghost>::~OctantCollection()
{
    p8est_destroy(forest);
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

} //end namespace

#endif