cmake_minimum_required(VERSION 3.18.0)

set(ALBRTLib_sources "" CACHE INTERNAL "")
message(STATUS ${ALBRTLib_sources})

add_subdirectory(lbm)
add_subdirectory(raja)
add_subdirectory(mpi)

message(STATUS ${ALBRTLib_sources})

blt_add_library(NAME ALBRTLib
                HEADERS ALBRT.h
                SOURCES "${ALBRTLib_sources}"
                INCLUDES "."
                DEPENDS_ON P4EST::P4EST SC::SC)

blt_print_target_properties(TARGET ALBRTLib)