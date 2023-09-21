#include <iostream>
#include <gtest/gtest.h>

#include "ALBRT.h"

using namespace ALBRT;

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    MPIManager::initialize();

    return RUN_ALL_TESTS();
    
    MPIManager::finalize();
}

TEST(OctantCollectionTestCase, OctantCollectionInitialization) {
    using T = double;
    using Lattice = D3Q19<T>;

    OctantCollection<T, Lattice, 20, 24, 28, 2> collection(4, 5, 6, 1.0, 3);
    OctantCollectionVTKManager<T, Lattice, 20, 24, 28, 2> vtkManager(collection);
    vtkManager.writeForest(0);
}
