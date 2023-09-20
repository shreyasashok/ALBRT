#ifndef OCTANT_COLLECTION_H
#define OCTANT_COLLECTION_H

namespace ALBRT {

/**
 * This class manages the forest of octants, manages device memory that
 * allows for operations on all process-local octants, and interfaces with
 * p4est for fringe exchange and adaption.
 * 
 * The template parameters are used to initialize the octants themselves.
 * 
 * T: data type of LBM storage quantities
 * Lattice: lattice descriptor
 * nx, ny, nz: dimensions of octant, including ghost layer.
 * nghost: number of ghost layers (must be even and at least 2)
*/
template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
class OctantCollection {

private:

    T octantPhysXLength; //physical x-dimension length of an octant in meters, not including ghost layers
    int maxLevel; //max level of refinement allowed

    


}; 

} //end namespace

#endif