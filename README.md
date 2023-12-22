# LarsieMPM

A General Finite Strain Elasto-Viscoplastic Material Point Method Framework in C++.

### Features

* **2D** and **3D**

* Symplectic Euler adaptive time stepping

* Adaptive background grid extension/reduction to particle domain

* Particle <-> grid interplation with **quadratic** and **cubic B-splines**

* Supports a linearly weighted combination of **PIC** and **FLIP**

* Supoorts **APIC** and **AFLIP**

* A nonlocal plasticity scheme is available for regularizing problems involving strain localization

* The following **elasticity models** are currently available:
    * Neo-Hookean model
    * Hencky's elasticity model


* The following **(visco)plastic models** are currently available and compatible with Hencky's elasticity model:
    * Von Mises
    * Drucker-Prager
    * Modified Cam-Clay
    * Perzyna-Von Mises
    * Perzyna-Drucker-Prager
    * Perzyna-Modified Cam-Clay
    * mu(I)-Drucker-Prager
    * mu(I)-Modified Cam-Clay
    * Implementing other (visco)plastic models is easy due to the general framework of the code


* Supports boundaries/objects formulated analytically as levelsets

* The supported boundary conditions are (the last two requiring a Coulomb friction parameter to be specified)
    1) **sticky**
    2) **slipping**
    3) **separating**


* Supports parallelization on shared memory with **OpenMP**

* Initial particle positions can be sampled using the Poisson disk sampling scheme by R. Bridson, ACM SIGGRAPH 2007, here based on the implementation by [Tommy Hinks](https://github.com/thinks/poisson-disk-sampling)

### Get started

1. Set up your simulation parameters and initial state in `mpm.cpp`   

2. Create a build directory: `mkdir build`

3. From the build directory (`cd build`), specify CMake options: `cmake -DCMAKE_BUILD_TYPE=Release ..`

4. Compile with `make -j <number of cores>`

5. Run the executable: `./src/mpm`

### Dependencies

The only required dependencies are **[CMake](https://cmake.org/)**, **[OpenMP](https://www.openmp.org/)** and the linear algebra library **[Eigen](https://eigen.tuxfamily.org/)**.

On Mac, you can install OpenMP through Homebrew with 

`brew install libomp`    

You might need to use the following CMake options depending on your version:   

`cmake -DCMAKE_BUILD_TYPE=Release -DOpenMP_CXX_FLAG="-Xclang -fopenmp" -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_CXX_LIB_NAMES=libomp -DOpenMP_C_FLAG="-Xclang -fopenmp" -DOpenMP_C_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_C_LIB_NAMES=libomp -DOpenMP_libomp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib ..`


### Output data

The location of the output data is specified by the user.
Particle data is saved as binary PLY-files (using [tinyply](https://github.com/ddiakopoulos/tinyply)) with the format (`out_part_X.ply`) where X represents the frame number (from 0 to `end_frame` as specified by the user).
This data format can be efficiently processed by SideFX's Houdini for visualization.
The optional outputting of the grid data saves only the positions, velocities and mass of the grid nodes in a CSV-file for each frame (`out_grid_X.csv`).
