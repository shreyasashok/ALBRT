#include <iostream>
#include <gtest/gtest.h>

#include "ALBRT.h"

using namespace ALBRT;

int main(int argc, char** argv){
    int retVal = 0;

    ::testing::InitGoogleTest(&argc, argv);
    ALBRTMain::initialize();

    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (MPIManager::rank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    int rankRetVal = RUN_ALL_TESTS();
    sc_MPI_Reduce(&rankRetVal, &retVal, 1, sc_MPI_INT, sc_MPI_MAX, 0, MPIManager::mpiCommWorld);
    
    ALBRTMain::finalize();
    return retVal;
}

TEST(OctantCollectionTestCase, OctantCollectionInitialization) {
    using T = double;
    using Lattice = D3Q19<T>;

    OctantCollection<T, Lattice, 20, 24, 28, 2> collection(3, 3, 3, 1.0, 3);
    OctantCollectionVTKManager<T, Lattice, 20, 24, 28, 2> vtkManager(collection);
    vtkManager.writeForest(0);
}
