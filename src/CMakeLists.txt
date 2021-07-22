if (USE_ABC)
	set(ALEMBIC_INCLUDES "/home/blatny/alembic/lib/")
	set(OPENEXR_INCLUDES "/usr/include/OpenEXR")
	include_directories("${ALEMBIC_INCLUDES}")
	include_directories("${OPENEXR_INCLUDES}")
	set(ALEMBIC_LIBRARIES "/home/blatny/alembic/build/lib/Alembic/libAlembic.so")
endif (USE_ABC)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if (USE_OMP)
	add_executable(mpm
				   mpm.cpp
			   	   simulation.cpp
				   tools.cpp
			       deformation_update.cpp
			       deformation_update_omp.cpp
				   plasticity.cpp
			       p2g.cpp
			       p2g_omp.cpp
			       g2p.cpp
			       g2p_omp.cpp
			       explicit_euler_update.cpp
				   explicit_euler_update_omp.cpp
				   p2g_nonlocal.cpp
				   g2p_nonlocal.cpp
				   plasticity_projection.cpp
				   remesh.cpp
			       )
	find_package(OpenMP REQUIRED)
	target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${ALEMBIC_LIBRARIES})
else (USE_OMP)
	add_executable(mpm
				   mpm.cpp
			   	   simulation.cpp
				   tools.cpp
			       deformation_update.cpp
				   plasticity.cpp
			       p2g.cpp
			       g2p.cpp
			       explicit_euler_update.cpp
				   p2g_nonlocal.cpp
				   g2p_nonlocal.cpp
				   plasticity_projection.cpp
				   remesh.cpp
			       )
	target_link_libraries (mpm Eigen3::Eigen ${ALEMBIC_LIBRARIES})
endif (USE_OMP)


# -lAlembic -L/home/blatny/alembic/build/lib/Alembic -I/home/blatny/alembic/lib/ -I/usr/include/OpenEXR
