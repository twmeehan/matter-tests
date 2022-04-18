cmake_minimum_required(VERSION 3.9)
SET (CMAKE_CXX_STANDARD 14)
project(Larsie)
add_subdirectory(src)

option (USE_ABC "Use Alembic part of the code" OFF)
option (USE_OMP "Use OpenMP part of the code" OFF)
