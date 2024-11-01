# Ma++er

Ma++er (pronounced _Matter_) is an abbreviation for *MPM in C++ for Elasto-viscoplastic **R**heologies

The Material Point Method (MPM) is a powerful numerical continuum method aimed at large deformation problems, typically attributed to Sulsky et al. (1994).

IMAGES HERE

## Features

* Both **2D** and **3D**


* Supports Update Stress Last **(USL)** and Modified Update Stress Last **(MUSL)** udated Lagrangian **symplectic** (explicit) MPM


* Particle <-> grid interplation with **quadratic** and **cubic B-splines** with various transfer schemes:
    *  Linearly weighted combination of **PIC** and **FLIP** (Stomakhin et al., 2013)
    *  **APIC** (Jiang et al., 2015)
    *  **AFLIP** (Fei et al., 2021)


* Adaptive background **grid extension/reduction** to particle domain


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


* Supports analytical objects and complex terrains


* The supported **boundary conditions** are (the last two requiring a **Coulomb friction parameter** to be specified)   
    * No-slip (in code called `STICKY`)
    * Separate slip (in code called `SEPARATE`)
    * Sticky slip (in code called `SLIP`, currently only available for plate objects)



* **Parallelized** on shared memory with [OpenMP](https://www.openmp.org/)

* Initial particle positions can be sampled using the **Poisson disk sampling** scheme by R. Bridson, ACM SIGGRAPH 2007, here based on the implementation by [Tommy Hinks](https://github.com/thinks/poisson-disk-sampling)

## Usage

#### How to run the code

1. Set up your simulation parameters and initial state in `mpm.cpp`. The default `mpm.cpp` file in the master branch sets up a simple granular collapse and explains the main options. In the `examples` folder, other examples are presented (to use one of these examples, simply copy it into the `src` folder and rename it `mpm.cpp`).

2. Create a build directory: `mkdir build`

3. From the build directory (`cd build`), specify CMake options: `cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VDB=ON ..`

4. Compile with `make -j <number of cores for compliation>` (NB: the number of threads for the _simulation_ is specified in `mpm.cpp`)

5. Run the executable: `./src/mpm`


#### Objects and terrains

Rigid objects and terrains (boundaries) are either   
    * formulated analytically as level sets (signed distance functions)   
    * or imported as `.vdb` level sets files using [OpenVDB](https://www.openvdb.org/)   

Analytical objects can be specified as a derived class from the general `ObjectGeneral` class. An example of this is `ObjectBump` which provides the terrain of a smooth bump used in the flow experiments in Viroulet et al. (2017). For the very common case of an axis-aligned plate, an `ObjectPlate` class has been made separate from `ObjectGeneral` class for speed and conveniance. In `ObjectPlate`, you can also assign a speed (and controls on the time-evokution of the speed) of the plate. Any plate must either a `top`, `bottom`, `front`, `back`, `left` or `right`. Its usage can be seen in the default example in `mpm.cpp` where it is used to simple set the ground (y = 0).

A terrain/object from a `.vdb` is stored in an instance of the `ObjectVdb` class which is also derived from `ObjectGeneral`. An example of a `.vdb` terrain of ... is found in the folder `levelsets`.

Note that all `ObjectGeneral` instances must be added to the std::vector `objects` and `ObjectPlate` instances are added to the std::vector `plates`.


#### Saving simulation data

The directory to save the output data is specified by the user in `mpm.cpp`.
Particle data is saved as **binary PLY-files** (using [tinyply](https://github.com/ddiakopoulos/tinyply)) with the format (`particles_fX.ply`) where X represents the frame number (from 0 to `end_frame` as specified by the user).

We recommend [SideFX's Houdini](https://www.sidefx.com) for visualization the particle data. In the file `visualize.hipnc`, we show how to make a simple visualization of the data. Some simple postprocesing can also be done directly in Houdini. PLY files can also be easily imported in Python for postprocessing. This is shown in the file `load_ply.py`.

Optionally, the grid data can also be saved if `save_grid = true`. Then, the grid data is saved as `grid_fX.ply`. By default, the grid data is not saved.

## Dependencies

The only required dependencies are **[CMake](https://cmake.org/)**, **[OpenMP](https://www.openmp.org/)** and the C++ linear algebra template library **[Eigen](https://eigen.tuxfamily.org/)**. The option `-DUSE_VDB=ON` also requires **[OpenVDB](https://www.openvdb.org/)**, however, this can be turned off if only analytic objects are used.

On Mac, you can install OpenMP through Homebrew with
`brew install libomp`
and Eigen can be obtained through
`brew install eigen`.

## Troubleshooting

* If OpenMP was installed through Homebrew, you might need to use the following CMake options when building on Mac, depending on your version:   
`cmake -DCMAKE_BUILD_TYPE=Release -DOpenMP_CXX_FLAG="-Xclang -fopenmp" -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_CXX_LIB_NAMES=libomp -DOpenMP_C_FLAG="-Xclang -fopenmp" -DOpenMP_C_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_C_LIB_NAMES=libomp -DOpenMP_libomp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib ..`


## Coming soon

* Option to export output files in other formats, notably `vtk`
* More efficient grid handling with [NanoVDB](https://www.openvdb.org/documentation/doxygen/NanoVDB_MainPage.html).


## License and attribution
Ma++er is licensed under ...
If you are interested in using Ma++er in commercial products or services, please do not hesitate to contact [Lars Blatny](https://larsblatny.github.io/) (lars.blatny@slf.ch).

If you use Ma++er in your research, please consider to cite works where this code has been used:   
Blatny, L., Gray, J.M.N.T. and Gaume, J. (2024) _A critical state μ(I)-rheology model for cohesive granular flows_, Journal of Fluid Mechanics, 997, p. A67. [doi:10.1017/jfm.2024.643](https://doi.org/10.1017/jfm.2024.643)

## Help?

Please contact [Lars Blatny](https://larsblatny.github.io/) (lars.blatny@slf.ch)
