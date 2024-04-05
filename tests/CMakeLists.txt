
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(
  elasticity
  elasticity.cpp
  ../src/tools.cpp
  ../src/simulation.cpp
  ../src/save_data.cpp
  ../src/remesh.cpp
  ../src/deformation_update_omp.cpp
  ../src/plasticity.cpp
  ../src/p2g_omp.cpp
  ../src/g2p_omp.cpp
  ../src/explicit_euler_update_omp.cpp
  ../src/boundary_collision.cpp
  ../src/pbc.cpp
)
target_link_libraries(
  elasticity
  GTest::gtest_main
  OpenMP::OpenMP_CXX
  Eigen3::Eigen
)

include(GoogleTest)
gtest_discover_tests(elasticity)
