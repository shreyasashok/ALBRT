cmake_minimum_required(VERSION 3.18.0)

set(OctantCollectionTest_depends ALBRTLib RAJA gtest)

if (ENABLE_CUDA)
    set(OctantCollectionTest_depends ${OctantCollectionTest_depends} cuda)
endif()

blt_add_executable(NAME OctantCollectionTest    SOURCES octantCollectionTest.cpp 
                                                DEPENDS_ON ${OctantCollectionTest_depends}
                                                OUTPUT_DIR "${CMAKE_BINARY_DIR}/bin/tests/OctantCollectionTest/"
                                                )
blt_add_test(NAME OctantCollectionTest COMMAND OctantCollectionTest)