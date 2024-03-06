cmake_minimum_required(VERSION 3.9)
SET (CMAKE_CXX_STANDARD 17)
project(Larsie)
add_subdirectory(src)

option (USE_VDB "Use OpenVDB part of the code" OFF)
