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
    * Peric Von Mises
    * Peric Drucker-Prager
    * Modified/Cohesive Cam Clay
    * Cohesive Quadratic
    * Implementing other (visco)plastic models is easy due to the general framework of the code


* Analytic objects formulated as levelsets are supported with 1) **sticky**, 2) **slipping** or 3) **separating** boundary conditions. Currently, only infinite plate objects are implemented

* Supports parallelization on shared memory with **OpenMP**

* A Python script that generates initial particle positions distributed according to the Poisson disk sampling scheme by R. Bridson, ACM SIGGRAPH 2007 is provided

### Get started

1. Set up your simulation parameters and initial state in `mpm.cpp`   

2. Create a build directory: `mkdir build`

3. From the build directory (`cd build`), specify CMake options: `cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=True ..`

4. Compile with `make -j <number of cores>`

5. Run the executable: `./src/mpm`

### Dependencies

The only required dependencies are **CMake** and the linear algebra library **Eigen**.

The parallelized version with **OpenMP** is optional. To use this version, make sure to set the CMake option `-DUSE_OMP=True`, and select the appropriate parallelized functions you want to use in `simulation.cpp`.

### Output data

The location of the output data is specified by the user. The data is by default saved as csv-files with the format (x, y, z, vx, vy, vz, ...) for both particles (`out_part_X.csv`) and grid (`out_grid_X.csv`) data where X represents the frame number (from 0 to `end_frame`).

### Validation

The code offers the possibility for a user-defined external force which may depend on the Lagrangian coordinates of the particles. This allows for the creation of manufactured solutions, which can be used to validate the code in the pure elastic case.
