cmake_minimum_required(VERSION 3.18.0)

list(APPEND ALBRTLib_sources "${CMAKE_CURRENT_SOURCE_DIR}/vtkMPIManager.cpp")
set(ALBRTLib_sources ${ALBRTLib_sources} PARENT_SCOPE)