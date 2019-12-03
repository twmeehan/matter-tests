add_executable(mpm mpm.cpp simulation.cpp)

set(ALEMBIC_INCLUDES "/home/blatny/alembic/lib/")
set(OPENEXR_INCLUDES "/usr/include/OpenEXR")
include_directories("${ALEMBIC_INCLUDES}")
include_directories("${OPENEXR_INCLUDES}")
set(ALEMBIC_LIBRARIES "/home/blatny/alembic/build/lib/Alembic/libAlembic.so")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (mpm Eigen3::Eigen ${ALEMBIC_LIBRARIES})


# -lAlembic -L/home/blatny/alembic/build/lib/Alembic -I/home/blatny/alembic/lib/ -I/usr/include/OpenEXR
