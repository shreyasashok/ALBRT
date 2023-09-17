#ifndef OCTANT_HPP
#define OCTANT_HPP

#include "octant.h"
#include "../raja/memoryManager.h"

namespace ALBRT {

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
Octant<T, Lattice, nx, ny, nz, nghost>::Octant()
{
    static_assert(nghost >= 2, "At least 2 ghost layers required!");
    static_assert(nghost % 2 == 0, "Number of ghost layers must be even!");


    // Allocate memory for this octant.
    size_t lbmBlockAllocSize =  (nx + 1) *  // This accounts for LBM populations.
                                (ny + 2) *  // The +1/+2 accounts for the extra space
                                (nz + 2);   // required by the EsoPull streaming scheme.

    size_t macroSize = nx*ny*nz;    // This accounts for macroscopic quantities.
                                    // (Density, Velocity, Force)
                                    // No need for the extra +1/+2 here.

    size_t blockDataSize =  Lattice::q * lbmBlockAllocSize +   
                            (1 + 2 * Lattice::d) * macroSize;    
                            

    rawData = MemoryManager::allocDevice<T>(blockDataSize); // Allocate the raw data

    evenData = MemoryManager::allocDeviceManaged<T*>(Lattice::dataSize); // Allocate the even and odd 
                                                                           // pointer-to-pointer structures.
    oddData = MemoryManager::allocDeviceManaged<T*>(Lattice::dataSize);

    // Handle the zeroeth population:
    evenData.get()[0] = rawData.get();
    oddData.get()[0] = rawData.get();

    // Loop over remaining populations and set up EsoPull pointer shifts.
    for (int iPop = 1; iPop <= Lattice::q/2; iPop++) {
        int dx = Lattice::c(iPop, 0);
        int dy = Lattice::c(iPop, 1);
        int dz = Lattice::c(iPop, 2);

        int shift = getCellIndex(dx, dy, dz);
        if (shift < 0) {
            evenData.get()[iPop] = rawData.get() + lbmBlockAllocSize*iPop - shift;
            oddData.get()[Lattice::opp(iPop)] = evenData.get()[iPop] + shift;

            evenData.get()[Lattice::opp(iPop)] = rawData.get() + lbmBlockAllocSize*Lattice::opp(iPop);
            oddData.get()[iPop] = evenData.get()[Lattice::opp(iPop)] - shift;
        }
        else {
            evenData.get()[iPop] = rawData.get() + lbmBlockAllocSize*iPop;
            oddData.get()[Lattice::opp(iPop)] = evenData.get()[iPop] + shift;

            evenData.get()[Lattice::opp(iPop)] = rawData.get() + lbmBlockAllocSize*Lattice::opp(iPop) + shift;
            oddData.get()[iPop] = evenData.get()[Lattice::opp(iPop)] - shift;
        }
    }

    // Loop over macroscopic variables and assign their pointer locations.
    for (int macro = Lattice::q; macro < Lattice::q + 1 + 2*Lattice::d; macro++) {
        evenData.get()[macro] = rawData.get() + lbmBlockAllocSize*Lattice::q + macroSize*macro;
        oddData.get()[macro]  = rawData.get() + lbmBlockAllocSize*Lattice::q + macroSize*macro;
    }

    // Initialize the cell mask (one integer per cell)
    cellMask = MemoryManager::allocDevice<int>(macroSize);
}

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
T* ALBRT::Octant<T, Lattice, nx, ny, nz, nghost>::getRawData()
{
    return rawData.get();
}

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
T** ALBRT::Octant<T, Lattice, nx, ny, nz, nghost>::getEvenData()
{
    return evenData.get();
}

template<typename T, class Lattice, int nx, int ny, int nz, int nghost>
T** ALBRT::Octant<T, Lattice, nx, ny, nz, nghost>::getOddData()
{
    return oddData.get();
}

}

#endif