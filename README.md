# LarsieMPM

An implementation of the Material Point Method (MPM) in the finite-strain elastoplastic framework.

### Features and limitations

Currently only supports 2D and quadratic splines.

Supports Neo-Hookean and St. Venant-Kirchhoff elasticity.

So far the only plastic model is Von Mises with Hencky strain, which can only be used with the St. Venant-Kirchhoff elastic model.

Implementing other elastic or plastic models is easy due to the general framework of the code.

### How to compile and run

`mkdir build`

`cd build`

`mkdir dumps`

`cmake -DCMAKE_BUILD_TYPE=Release ..`

Set up your simulation parameters and initial state in the `mpm.cpp` file.   

Compile from the `build` directory with `make`.  

Run from the `build` directory with the following command:  

`./src/mpm`

### Output data

The output is saved in the `dumps` directory as csv-files in the format (x, y, z, vx, vy, vz) for both particles (`out_part_X.csv`) and grid (`out_grid_X.csv`) data where X represents the frame number (from 0 to `end_frame`). The particle data file may contain more information than just position and velocity, e.g., pressure, plastic deviatoric/volumetric strain.

### Validation

The code offers the possibility for a user-defined external force which may depend on the Lagrangian coordinates of the particles. This allows for the creation of manufactured solutions, which can be used to validate the code (at least in the pure elastic case).

### Performance

Larsie is continously implemented according to the principle that premature optimization is the root of all evil. Basic optimizations (precomputations, clever particle-grid loops, memory reads, etc...) are currently still being explored before moving onto parallelization.

### Dependencies

Larsie only relies on the linear algebra library Eigen. 
