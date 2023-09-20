cmake_minimum_required(VERSION 3.18.0)


set(EsoPullTest_depends ALBRTLib RAJA gtest)

if (ENABLE_CUDA)
    set(EsoPullTest_depends ${EsoPullTest_depends} cuda)
endif()

message(STATUS "Here! ${EsoPullTest_depends}")
blt_add_executable(NAME EsoPullTest SOURCES esoPullTest.cpp DEPENDS_ON ${EsoPullTest_depends})

blt_print_target_properties(TARGET EsoPullTest)

blt_add_test(NAME EsoPullTest COMMAND EsoPullTest)

# target_include_directories(asdf PUBLIC "${CMAKE_SOURCE_DIR}/src" ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
# target_compile_definitions(asdf PUBLIC $ALBRT_DEFINITIONS)

