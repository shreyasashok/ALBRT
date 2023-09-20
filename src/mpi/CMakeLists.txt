cmake_minimum_required(VERSION 3.18.0)

list(APPEND ALBRTLib_sources "${CMAKE_CURRENT_SOURCE_DIR}/mpiManager.cpp")
message(STATUS ${ALBRTLib_sources})
set(ALBRTLib_sources ${ALBRTLib_sources} PARENT_SCOPE)