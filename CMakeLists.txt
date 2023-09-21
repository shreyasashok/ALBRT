cmake_minimum_required(VERSION 3.18.0)
project(ALBRT VERSION 0.1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 17)

if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++17" CACHE STRING "")
endif()

if (ENABLE_CUDA)
    set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda" CACHE PATH "")
    set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc" CACHE PATH "")
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE PATH "")
    set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    
    set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
    list(APPEND ALBRT_DEFINITIONS "ALBRT_USE_CUDA")
endif ()

include(blt/SetupBLT.cmake)

find_package(VTK COMPONENTS 
            CommonCore
            CommonDataModel
            IOLegacy
            IOXML
            IOParallelXML
            ParallelMPI
            )

find_package(SC CONFIG)
if(NOT SC_FOUND)
    find_package(SC REQUIRED)
endif()

find_package(P4EST CONFIG)
if(NOT P4EST_FOUND)
    find_package(P4EST REQUIRED)
endif()

find_package(RAJA REQUIRED)
# find_package(GTest CONFIG REQUIRED)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(ALBRTExec main.cpp)

