find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if (USE_OMP)
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
				   plasticity_projection.cpp
				   p2g_nonlocal.cpp
				   g2p_nonlocal.cpp
				   pbc.cpp
			       )
	find_package(OpenMP REQUIRED)
	target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${ALEMBIC_LIBRARIES})
else (USE_OMP)
	add_executable(mpm
				   mpm.cpp
			   	   simulation.cpp
				   save_data.cpp
			       deformation_update.cpp
				   plasticity.cpp
			       p2g.cpp
			       g2p.cpp
			       explicit_euler_update.cpp
				   plasticity_projection.cpp
				   remesh.cpp
				   validate_rma.cpp
				   boundary_collision.cpp
				   pbc.cpp
			       )
	target_link_libraries (mpm Eigen3::Eigen ${ALEMBIC_LIBRARIES})
endif (USE_OMP)
