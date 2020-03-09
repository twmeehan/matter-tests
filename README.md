# LarsieMPM

A simple implementation of the Material Point Method (MPM) in the elastoplastic framework.

### Features and limitations:

Currently only supports 2D and quadratic splines. Supports Neo-Hookean and St. Venant-Kirchhoff elasticity. So far the only plastic model is Von Mises with Hencky strain, which can only be used with the St. Venant-Kirchhoff elastic model.

### How to compile and run:

`mkdir build`

`cd build`

`mkdir dumps`

`cmake -DCMAKE_BUILD_TYPE=Release ..`

`make`

Run from the `build` folder with the following command:

`./src/mpm`

### Output data:

The output is saved in the `dumps` directory as csv-files in the format (x, y, z, vx, vy, vz) for both particles (`out_part_X.csv`) and grid (`out_grid_X.csv`) data where X represents the frame number (from 0 to `end_frame`).
