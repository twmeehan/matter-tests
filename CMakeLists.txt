cmake_minimum_required(VERSION 3.9)
SET (CMAKE_CXX_STANDARD 17)
project(matter)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

option (USE_VDB "Use OpenVDB part of the code" OFF)
