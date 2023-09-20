#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#include <sc.h>

namespace ALBRT {

struct MPIManager {

    static int rank;
    static int numRanks;

    static sc_MPI_Comm mpiCommWorld;

    /**
     * Initialize MPI through the SC library.
    */
    static void initialize();

    /**
     * Finalize MPI through the SC library.
    */
    static void finalize();
};

} //end namespace

#endif