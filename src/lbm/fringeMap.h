#ifndef FRINGE_MAPS_H
#define FRINGE_MAPS_H

#include "neighborRelationship.h"
#include "commLattice.h"
#include "RAJA/RAJA.hpp"

namespace ALBRT {

/**
 * Contains mapping functions which determine size and index mapping
 * for fringe exchanges, including match, coarse-to-fine, fine-to-coarse,
 * and boundary relationship types.
 * 
 * nghostBlock refers to how many ghost cells each structured block contains.
 * But we don't always have to transfer information to all the ghost cells.
 * To save communication time, we might only transfer a limited number of ghost
 * layers to the fringe.
 * nghostTransfer therefore controls how many ghost cells are actually transferred.
 * 
 * Typically, there will be a FringeMap for the populations and a separate FringeMap
 * for the macroscopic variables. For many boundary conditions, ghost data is only needed
 * for the macroscopic variables and not the populations themselves. By omitting layers 
 * of population ghost points, we can save greatly on communication.
*/
template <typename T, class Lattice, int nx, int ny, int nz, 
          int nghostBlock, int nghostTransfer>
struct FringeMap {

    /**
     * Get the number of elements required in the pack buffer to pack a face.
    */
    RAJA_DEVICE static int getFacePackSize(NeighborType type, int sendDirection, int varStart, int varEnd);

    /**
     * Pack a single cell into the pack buffer. This function is intended
     * to be called as part of a RAJA workgroup. 
    */
    RAJA_DEVICE static void packFace(NeighborType type, int sendDirection, int varStart, int varEnd,
                              T ** block, T* buffer, int * mask, int index);

    /**
     * Unpack a single cell from the pack buffer. This function is intended
     * to be called in parallel as part of a RAJA workgroup. 
    */
    RAJA_DEVICE void unpackFace(NeighborType type, int recvDirection, int varStart, int varEnd,
                              T ** block, T* buffer, int * mask, int index);

private:
    RAJA_DEVICE static int indexFromOriginAndSize(int regionIndex, int iXStart, int iYStart, int iZStart, int iXSize, int iYSize);

    static constexpr int nxi = nx - 2*nghostBlock; //nx of non-ghost region
    static constexpr int nyi = ny - 2*nghostBlock; //nx of non-ghost region
    static constexpr int nzi = nz - 2*nghostBlock; //nx of non-ghost region

    // First three numbers are origin, second three are size
    static constexpr faceFringeMapMatchPack[FaceCommLattice::q][6] = 
        {   {nghostBlock, nghostBlock, nghostBlock, nghostTransfer, nyi, nzi},
            {nx-nghostBlock-nGhostTransfer, nghostBlock, nghostBlock, nghostTransfer, nyi, nzi},
            {nghostBlock, nghostBlock, nghostBlock, nxi, nghostTransfer, nzi},
            {nghostBlock, ny-nghostBlock-nghostTransfer, nghostBlock, nxi, nghostTransfer, nzi},
            {nghostBlock, nghostBlock, nghostBlock, nxi, nyi, nghostTransfer},
            {nghostBlock, nghostBlock, nz-nghostBlock-nGhostTransfer, nxi, nyi, nghostTransfer}
        };

    static constexpr faceFringeMapMatchUnpack[FaceCommLattice::q][6] = 
        {   {nghostBlock-nghostTransfer, nghostBlock, nghostBlock, nghostTransfer, nyi, nzi},
            {nx-nghostBlock, nghostBlock, nghostBlock, nghostTransfer, nyi, nzi},
            {nghostBlock, nghostBlock-nGhostTransfer, nghostBlock, nxi, nghostTransfer, nzi},
            {nghostBlock, ny-nghostBlock, nghostBlock, nxi, nghostTransfer, nzi},
            {nghostBlock, nghostBlock, nghostBlock-nghostTransfer, nxi, nyi, nghostTransfer},
            {nghostBlock, nghostBlock, nz-nghostBlock, nxi, nyi, nghostTransfer}
        };
};

} //end namespace

#endif