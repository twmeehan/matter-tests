find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(mpm
			   mpm.cpp
			   tools.cpp
			   plastic_models.cpp
		   	   simulation.cpp
			   save_data.cpp
			   remesh.cpp
		       deformation_update_omp.cpp
			   plasticity.cpp
		       p2g_omp.cpp
		       g2p_omp.cpp
			   explicit_euler_update_omp.cpp
			   boundary_collision.cpp
			   validate_rma.cpp
			   pbc.cpp
		       )
find_package(OpenMP REQUIRED)
target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${ALEMBIC_LIBRARIES})
