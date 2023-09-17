cmake_minimum_required(VERSION 3.0.0)

add_subdirectory(lbm)
add_subdirectory(raja)

blt_add_library(NAME ALBRTLib
                HEADERS ALBRT.h
                INCLUDES ".")

blt_print_target_properties(TARGET ALBRTLib)

message(STATUS "Done this!")