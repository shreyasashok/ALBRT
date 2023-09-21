#include "mpiManager.h"

namespace ALBRT {

sc_MPI_Comm MPIManager::mpiCommWorld;
int MPIManager::rank;
int MPIManager::numRanks;

void MPIManager::initialize()
{
    int mpiret = sc_MPI_Init(NULL, NULL);
    SC_CHECK_MPI(mpiret);
    mpiCommWorld = sc_MPI_COMM_WORLD;

    // init SC with production level logging
    sc_init(mpiCommWorld, 1, 1, NULL, SC_LP_PRODUCTION);

    sc_MPI_Comm_rank(mpiCommWorld, &rank);
    sc_MPI_Comm_size(mpiCommWorld, &numRanks);
}

void MPIManager::finalize()
{
    sc_finalize();
    int mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);
}

}