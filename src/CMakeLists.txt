cmake_minimum_required(VERSION 3.18.0)

set(ALBRTLib_sources "" CACHE INTERNAL "")

add_subdirectory(lbm)
add_subdirectory(raja)
add_subdirectory(mpi)
add_subdirectory(vtk)

list(APPEND ALBRTLib_sources "${CMAKE_CURRENT_SOURCE_DIR}/ALBRT.cpp")

blt_add_library(NAME ALBRTLib
                HEADERS ALBRT.h
                SOURCES "${ALBRTLib_sources}"
                INCLUDES "."
                DEPENDS_ON  P4EST::P4EST SC::SC RAJA
                            VTK::CommonCore
                            VTK::CommonDataModel
                            VTK::IOLegacy
                            VTK::IOXML
                            VTK::IOParallelXML
                            VTK::ParallelMPI )

blt_print_target_properties(TARGET ALBRTLib)