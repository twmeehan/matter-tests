# Larsie

An implementation of the Material Point Method (MPM) in the finite-strain elastoplastic framework.

### Features and limitations

Supports **2D** and **3D**.

Supports **quadratic and cubic B-splines**, and **PIC**, **FLIP** and **PIC-FLIP** for particle-grid interpolation.

Supports **Neo-Hookean** and **St. Venant-Kirchhoff** elasticity.

Supports **Von Mises** and **Drucker-Prager** plasticity, which are to be used with the St. Venant-Kirchhoff elastic model. Both the Von Mises model and Drucker-Prager model has an optional **strain-softening** feature.

Implementing other elastic or plastic models is easy due to the general framework of the code.

**Analytic objects** (formulated as levelsets) are supported with **sticky** or **slipping** (with user-defined **friction**) boundary conditions. Currently, only infinite plate objects are implemented.

### How to compile and run

`mkdir build`

`cd build`

`cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=True ..`

Set up your simulation parameters and initial state in the `mpm.cpp` file.   

Compile from the `build` directory with `make`.  

Run from the `build` directory with the following command:  

`./src/mpm`

### Output data

Assuming the code is run from the `build` directory, the output is saved in the directory `dumps/<sim_name>` where `sim_name` is specified in the setup. The data is saved as csv-files with the format (x, y, z, vx, vy, vz, ...) for both particles (`out_part_X.csv`) and grid (`out_grid_X.csv`) data where X represents the frame number (from 0 to `end_frame`).

### Pre- and post processing
Two python scripts are provided for preprocessing and postprocessing, respectively. The preprocessing script generates samples distributed according to the Poisson Disk Sampling by R. Bridson, ACM SIGGRAPH 2007.

### Validation

The code offers the possibility for a user-defined external force which may depend on the Lagrangian coordinates of the particles. This allows for the creation of manufactured solutions, which can be used to validate the code (at least in the pure elastic case). The plastic models can easily be validated by plotting *p* (pressure) and *q* (Mises equivalent stress) for each particle in time.

### Performance

Larsie is continously implemented according to the principle *premature optimization is the root of all evil*. Basic optimizations (precomputations, clever particle-grid loops, memory reads, etc...) are still being explored. Although still under development, Larsie also supports **parallelization** on shared memory with **OpenMP**.

### Dependencies

Larsie only relies only on CMake and the linear algebra library Eigen.
To use OpenMP, use the CMake option `-DUSE_OMP=True`, and make sure to select the appropriate parallelized functions you want to use in `simulation.cpp`.
