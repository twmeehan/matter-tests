# LarsieMPM

A General Finite Strain Elasto-Viscoplastic Material Point Method Framework in C++.

### Features

* **2D** and **3D**

* Explicit Euler adaptive time stepping

* Adaptive background grid extension/reduction to particle domain

* Supports a linearly weighted combination of **PIC** and **FLIP** with **quadratic** and **cubic B-splines** for particle-grid interplation

* A nonlocal plasticity scheme is available for regularizing problems involving strain localization


* The following **elasticity models** are currently available:
    * Neo-Hookean
    * Saint Venant-Kirchhoff with Hencky strain


* The following **(visco)plastic models** are currently available and compatible with the SVK-Hencky model above:
    * Von Mises
    * Drucker-Prager
    * Modified Cam-Clay
    * Perzyna-Von Mises
    * Perzyna-Drucker-Prager
    * Perzyna-Modified Cam-Clay
    * mu(I)-Drucker-Prager
    * mu(I)-Modified Cam-Clay
    * Implementing other (visco)plastic models is easy due to the general framework of the code


* Analytic objects formulated as levelsets are supported with 1) **sticky**, 2) **slipping** or 3) **separating** boundary conditions. Currently, the only implemented objects are axis-aligned plates which can move with a prescribed velocity.

* Supports parallelization on shared memory with **OpenMP**

* Initial particle positions can be sampled using the Poisson disk sampling scheme by R. Bridson, ACM SIGGRAPH 2007, here based on the implementation by [Tommy Hinks](https://github.com/thinks/poisson-disk-sampling)

### Get started

1. Set up your simulation parameters and initial state in `mpm.cpp`   

2. Create a build directory: `mkdir build`

3. From the build directory (`cd build`), specify CMake options: `cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=True ..`

4. Compile with `make -j <number of cores>`

5. Run the executable: `./src/mpm`

### Dependencies

The only required dependencies are **[CMake](https://cmake.org/)** and the linear algebra library **[Eigen](https://eigen.tuxfamily.org/)**.

The parallelized version with **[OpenMP](https://www.openmp.org/)** is optional. To use this version, make sure to set the CMake option `-DUSE_OMP=True`, and select the appropriate parallelized functions you want to use in `simulation.cpp`.

### Output data

The location of the output data is specified by the user.
Particle data is saved as binary PLY-files (using [tinyply](https://github.com/ddiakopoulos/tinyply)) with the format (`out_part_X.ply`) where X represents the frame number (from 0 to `end_frame` as specified by the user).
This data format can be efficiently processed by SideFX's Houdini for visualization.
The optional outputting of the grid data saves only the positions, velocities and mass of the grid nodes in a CSV-file for each frame (`out_grid_X.csv`).


### Validation

The code offers the possibility for a user-defined external force which may depend on the Lagrangian coordinates of the particles. This allows for the creation of manufactured solutions, which can be used to validate the code in the pure elastic case.
