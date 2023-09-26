#ifndef FRINGE_MAPS_HH
#define FRINGE_MAPS_HH

#include "fringeMap.h"

namespace ALBRT {

template<typename T, class Lattice, int nx, int ny, int nz, int nghostBlock, int nghostTransfer>
int FringeMap<T, Lattice, nx, ny, nz, nghostBlock, nghostTransfer>::
getFacePackSize(NeighborType type, int sendDirection, int varStart, int varEnd)
{
    if (type == MATCH) {
        return  faceFringeMapMatchPack[sendDirection][3] *
                faceFringeMapMatchPack[sendDirection][4] *
                faceFringeMapMatchPack[sendDirection][5] *
                (varEnd-varStart);
    }
    else if (type == FINE) {
        return -1; //TODO
    }
    else if (type == COARSE) {
        return -1; //TODO
    }
}

template<typename T, class Lattice, int nx, int ny, int nz, int nghostBlock, int nghostTransfer>
void FringeMap<T, Lattice, nx, ny, nz, nghostBlock, nghostTransfer>::
packFace(NeighborType type, int sendDirection, int varStart, int varEnd,
         T ** block, T* buffer, int * mask, int index)
{
    const T nanVal = sqrt(T(-1.0));

    if (type == MATCH) {
        int sourceIndex = indexFromOriginAndSize(index,
                            faceFringeMapMatchPack[sendDirection, 0],
                            faceFringeMapMatchPack[sendDirection, 1],
                            faceFringeMapMatchPack[sendDirection, 2],
                            faceFringeMapMatchPack[sendDirection, 3],
                            faceFringeMapMatchPack[sendDirection, 4]);

        int oneVarPackSize =    faceFringeMapMatchPack[sendDirection][3] *
                                faceFringeMapMatchPack[sendDirection][4] *
                                faceFringeMapMatchPack[sendDirection][5];

        
        for (int iVar = varStart; iVar <= varEnd; iVar++) {
            buffer[index + oneVarPackSize*(iVar-varStart)] = 
                block[iVar][sourceIndex];
        }
    }
}

template<typename T, class Lattice, int nx, int ny, int nz, int nghostBlock, int nghostTransfer>
void FringeMap<T, Lattice, nx, ny, nz, nghostBlock, nghostTransfer>::
unpackFace(NeighborType type, int recvDirection, int varStart, int varEnd,
         T ** block, T* buffer, int * mask, int index)
{
    const T nanVal = sqrt(T(-1.0));

    if (type == MATCH) {
        int destIndex = indexFromOriginAndSize(index,
                            faceFringeMapMatchUnpack[recvDirection, 0],
                            faceFringeMapMatchUnpack[recvDirection, 1],
                            faceFringeMapMatchUnpack[recvDirection, 2],
                            faceFringeMapMatchUnpack[recvDirection, 3],
                            faceFringeMapMatchUnpack[recvDirection, 4]);

        int oneVarPackSize =    faceFringeMapMatchUnpack[recvDirection][3] *
                                faceFringeMapMatchUnpack[recvDirection][4] *
                                faceFringeMapMatchUnpack[recvDirection][5];

        
        for (int iVar = varStart; iVar <= varEnd; iVar++) {
            block[iVar][destIndex] = 
                buffer[index + oneVarPackSize*(iVar-varStart)] 
        }
    }
}

template <typename T, class Lattice, int nx, int ny, int nz, int nghostBlock, int nghostTransfer>
int FringeMap<T, Lattice, nx, ny, nz, nghostBlock, nghostTransfer>::
indexFromOriginAndSize(int regionIndex, int iXStart, int iYStart, int iZStart, int iXSize, int iYSize)
{ 
    int indices[3];
    indices[0] = regionIndex % iXSize + iXStart;
    indices[1] = regionIndex/iXSize % iYSize + iYStart;
    indices[2] = regionIndex/(iXSize*iYSize) + iZStart;

    return indices[0] + indices[1]*nx + indices[2]*nx*ny;
}

} //end namespace

#endif