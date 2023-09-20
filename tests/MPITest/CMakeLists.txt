cmake_minimum_required(VERSION 3.18.0)

set(MPITest_depends ALBRTLib RAJA gtest)

if (ENABLE_CUDA)
    set(MPITest_depends ${MPITest_depends} cuda)
endif()

blt_add_executable(NAME MPITest SOURCES mpiTest.cpp DEPENDS_ON ${MPITest_depends})
blt_add_test(NAME MPITest COMMAND MPITest)