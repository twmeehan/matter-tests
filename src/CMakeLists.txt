find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(mpm
			   mpm.cpp
			   tools.cpp
			   simulation.cpp
			   save_data.cpp
			   remesh.cpp
			   deformation_update.cpp
			   plasticity.cpp
			   p2g.cpp
			   g2p.cpp
			   explicit_euler_update.cpp
			   boundary_collision.cpp
			   update_dt.cpp
			   position_update.cpp
			   musl.cpp
			   pbc.cpp
			   )

if (USE_VDB)

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

	find_package(OpenVDB 3.1 QUIET)
	include_directories(${OPENVDB_INCLUDE_DIRS} SYSTEM)

	find_package(TBB REQUIRED)
	include_directories(${TBB_INCLUDE_DIRS} SYSTEM)

	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS} SYSTEM)

	if (APPLE)
		message("-- MESSAGE: Will not include Half libraries")
		target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${OPENVDB_LIBRARIES} ${OpenEXR_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES})
	else (APPLE)
		message("-- MESSAGE: Will include Half libraries")
		find_package(Half QUIET)
		include_directories(${HALF_INCLUDE_DIRS} SYSTEM)

		target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${OPENVDB_LIBRARIES} ${OpenEXR_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} ${HALF_LIBRARIES})
	endif (APPLE)

else (USE_VDB)

	target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen)

endif (USE_VDB)
