#ifndef OCTANT_COLLECTION_H
#define OCTANT_COLLECTION_H

#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_iterate.h>
#include <vector>
namespace ALBRT {

/**
 * Contains internal information about each octant.
 * This is used as "user data" in p4est. Some of the 
 * information is duplicated from the p4est information,
 * but this is for our convenience since they use different
 * conventions than we want to use.
*/
template <typename T>
struct OctantInternalInfo {
    int proc;
    int index;
    int level;

    T origin[3];
    T physLength[3];
};


/**
 * Forward declaration; VTK manager class is a friend, so it
 * can manipulate private members in OctantCollection.
*/
template <typename T, class Lattice, int nx, int ny, int nz, int nghost>
class OctantCollectionVTKManager;

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
    friend class OctantCollectionVTKManager<T, Lattice, nx, ny, nz, nghost>;

private:
    typedef typename Lattice::Real S;

    T baseOctantPhysLength[3]; //physical dimension length of an octant in meters, not including ghost layers
    int maxLevel; //max level of refinement allowed

    p8est_connectivity_t * connectivity;
    p8est_t              * forest;

    std::vector<OctantInternalInfo<T>> localOctantInfos;

    static int numLocalOctants;
    static void initOctantIndices(p8est_iter_volume_info_t * info, void *user_data);

public:
    /** 
     * Initializes the forest. We use the p4est "brick" connectivity for all grids.
     * The octant physical length corresponds to the X dimension. The other physical lengths
     * are determined based on the given octant dimensions.
     * 
     * The template parameters correspond to those required to initialize an octant (see octant.h).
     * 
     * Periodicity must be specified in the constructor, because this affects the p4est connectivity itself.
    */
    OctantCollection(int numOctantsX, int numOctantsY, int numOctantsZ, T octantPhysXLength, int maxLevel,
                     bool periodicX = false, bool periodicY = false, bool periodicZ = false);

    ~OctantCollection();

    

}; 

} //end namespace

#endif