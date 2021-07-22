#include "simulation.hpp"

void Simulation::remesh(){

    // ACTUAL min and max position of particles
    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = (*min_x_it)(0);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);


    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;
    T Lz = max_z - min_z;

    // safety_factor = 2 means we have a grid which has a grid point 2*dx from the boundary particle
    // Assuming local approach, the grid point 2dx away will not be be given a nonzero value
    unsigned int safety_factor = std::max((unsigned int)2, nonlocal_support);

    // NEW number of grid-dx's per side length
    Nx = std::ceil(Lx * one_over_dx) + 1 + 2*safety_factor;
    Ny = std::ceil(Ly * one_over_dx) + 1 + 2*safety_factor;
    Nz = std::ceil(Lz * one_over_dx) + 1 + 2*safety_factor;

    T low_x  = min_x - dx * safety_factor;
    T high_x = max_x + dx * safety_factor;
    T low_y  = min_y - dx * safety_factor;
    T high_y = max_y + dx * safety_factor;
    T low_z  = min_z - dx * safety_factor;
    T high_z = max_z + dx * safety_factor;

    debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x, high_x, Nx);
    grid.y = linspace(low_y, high_y, Ny);
    grid.z = linspace(low_z, high_z, Nz);

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
    grid.zc = grid.z[0];

    grid.v.resize(Nx*Ny*Nz);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(Nx*Ny*Nz); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(Nx*Ny*Nz); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}

void Simulation::remeshFixedInit(){

    // ACTUAL min and max position of particles
    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = (*min_x_it)(0);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);

    // Save for remeshFixedCont
    max_y_init = max_y;

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;
    T Lz = max_z - min_z;

    // safety_factor = 2 means we have a grid which has a grid point 2*dx from the boundary particle
    // Assuming local approach, the grid point 2dx away will not be be given a nonzero value
    unsigned int safety_factor = std::max((unsigned int)2, nonlocal_support);

    Nx = std::ceil(Lx * one_over_dx) + 1 + 2*safety_factor;
    Ny = std::ceil(Ly * one_over_dx) + 1 + 2*safety_factor;
    Nz = std::ceil(Lz * one_over_dx) + 1 + 2*safety_factor;

    // save for remeshFixedCont
    Ny_init = Ny;

    T low_x_init  = min_x - dx * safety_factor;
    T high_x_init = max_x + dx * safety_factor;
    low_y_init    = min_y - dx * safety_factor;
    high_y_init   = max_y + dx * safety_factor;
    T low_z_init  = min_z - dx * safety_factor;
    T high_z_init = max_z + dx * safety_factor;

    debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x_init, high_x_init, Nx);
    grid.y = linspace(low_y_init, high_y_init, Ny);
    grid.z = linspace(low_z_init, high_z_init, Nz);

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
    grid.zc = grid.z[0];

    grid.v.resize(Nx*Ny*Nz);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(Nx*Ny*Nz); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(Nx*Ny*Nz); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}

void Simulation::remeshFixedCont(){


    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);


    unsigned int reduction_factor = std::floor( std::max(0.0,(max_y_init-max_y)/dx - 1e-8*dx) );
    Ny       = Ny_init     - reduction_factor;
    T high_y = high_y_init - reduction_factor * dx;

    debug("               reduction factor =  ", reduction_factor);
    debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.y = linspace(low_y_init, high_y, Ny);

    grid.v.resize(Nx*Ny*Nz);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(Nx*Ny*Nz); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(Nx*Ny*Nz); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}
