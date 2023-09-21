#ifndef OCTANT_COLLECTION_VTK_MANAGER_H
#define OCTANT_COLLECTION_VTK_MANAGER_H

#include "../lbm/octantCollection.h"

#include <string>

namespace ALBRT {

/**
 * A VTK output manager for the octant collection.
 * This class will handle output of VTK files for any
 * kind of data related directly to the octant collection.
*/
template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
class OctantCollectionVTKManager {
    
    std::string outputDirectory;
    OctantCollection<T, Lattice, nx, ny, nz, nghost>& octCollection;

public:

    /**
     * Constructor; specify alternate output directory if preferred.
    */
    OctantCollectionVTKManager(OctantCollection<T, Lattice, nx, ny, nz, nghost>& octCollectionIn,
                                    std::string outputDirectory = "./vtk_out/");

    /**
     * Write a lightweight representation of the octree forest to VTK.
     * This function uses the .pvtu format.
    */
    void writeForest(long timestep, std::string filename = "forest");



};



} //end namespace

#endif