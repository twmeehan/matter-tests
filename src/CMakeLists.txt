add_executable(mpm mpm.cpp simulation.cpp)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (mpm Eigen3::Eigen)
