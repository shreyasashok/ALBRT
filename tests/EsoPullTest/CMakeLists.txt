cmake_minimum_required(VERSION 3.18.0)


set(EsoPullTest_depends ALBRTLib RAJA gtest)

if (ENABLE_CUDA)
    set(EsoPullTest_depends ${EsoPullTest_depends} cuda)
endif()

blt_add_executable(NAME EsoPullTest     SOURCES esoPullTest.cpp 
                                        DEPENDS_ON ${EsoPullTest_depends}
                                        OUTPUT_DIR "${CMAKE_BINARY_DIR}/bin/tests/EsoPullTest/")

blt_print_target_properties(TARGET EsoPullTest)

blt_add_test(NAME EsoPullTest COMMAND EsoPullTest)

