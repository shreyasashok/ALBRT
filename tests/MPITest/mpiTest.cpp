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

TEST(MPITestCase, MPITestCase) {
    SUCCEED();    
}
