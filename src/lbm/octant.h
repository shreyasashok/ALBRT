#ifndef OCTANT_H
#define OCTANT_H

#include "../raja/memoryManager.h"

namespace ALBRT {

/**
 * A single LBM octant in the forest.
 * 
 * T: data type of LBM storage quantities
 * Lattice: lattice descriptor
 * nx, ny, nz: dimensions of octant, including ghost layer.
 * nghost: number of ghost layers (must be even and at least 2)
 * 
 * This class sets up LBM memory for this octant, and implements EsoPull streaming scheme. 
 * 
*/
template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
class Octant {

private:
    std::shared_ptr<T> rawData = nullptr;
    std::shared_ptr<T*> evenData = nullptr;
    std::shared_ptr<T*> oddData = nullptr;

    std::shared_ptr<int> cellMask = nullptr;

public:
    /**
     * Construct an octant, allocate memory, and set up streaming scheme.
    */
    Octant();

    /**
     * Get raw octant data pointer.
    */
    T* getRawData();

    /**
     * Get pointer-to pointer for even timestep populations.
    */
    T** getEvenData();

    /**
     * Get pointer-to pointer for odd timestep populations.
    */
    T** getOddData();

    /**
     * Gets flattened cell index from Cartesian indices.
    */
    static inline int getCellIndex(int iX, int iY, int iZ) {
        return iX + iY*nx + iZ*nx*ny;
    }

    /**
     * Gets Cartesian indices from flattened cell index.
    */
    static inline void getCellIndices(int index, int indices[3]) {
        indices[0] = index % nx;
        indices[1] = index/nx % ny;
        indices[2] = index/(nx*ny);
    }
};

} //end namespace

#include "octant.hpp" //include source because this is a templated class
#endif
