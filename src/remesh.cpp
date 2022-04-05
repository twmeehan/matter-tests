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
#ifdef THREEDIM
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
#endif
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
#ifdef THREEDIM
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);
#endif

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;
#ifdef THREEDIM
    T Lz = max_z - min_z;
#endif

    // safety_factor = 2 means we have a grid which has a grid point 2*dx from the boundary particle
    // Assuming local approach, the grid point 2dx away will not be be given a nonzero value
    unsigned int safety_factor = std::max((unsigned int)2, nonlocal_support);

    // NEW number of grid-dx's per side length
    Nx = std::ceil(Lx * one_over_dx) + 1 + 2*safety_factor;
    Ny = std::ceil(Ly * one_over_dx) + 1 + 2*safety_factor;
#ifdef THREEDIM
    Nz = std::ceil(Lz * one_over_dx) + 1 + 2*safety_factor;
    grid_nodes = Nx*Ny*Nz;
#else
    grid_nodes = Nx*Ny;
#endif

    T low_x  = min_x - dx * safety_factor;
    T high_x = max_x + dx * safety_factor;
    T low_y  = min_y - dx * safety_factor;
    T high_y = max_y + dx * safety_factor;
#ifdef THREEDIM
    T low_z  = min_z - dx * safety_factor;
    T high_z = max_z + dx * safety_factor;
#endif

#ifdef WARNINGS
    #ifdef THREEDIM
        debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );
    #else
        debug("               grid   = (", Nx, ", ", Ny, ")"  );
    #endif
#endif

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x, high_x, Nx);
    grid.y = linspace(low_y, high_y, Ny);
#ifdef THREEDIM
    grid.z = linspace(low_z, high_z, Nz);
#endif

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
#ifdef THREEDIM
    grid.zc = grid.z[0];
#endif

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}

void Simulation::remeshFixedInit(unsigned int sfx, unsigned int sfy, unsigned int sfz){

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
#ifdef THREEDIM
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
#endif
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
#ifdef THREEDIM
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);
#endif

    // Save for remeshFixedCont
    max_x_init = max_x;
    min_x_init = min_x;

    max_y_init = max_y;
    min_y_init = min_y;

#ifdef THREEDIM
    max_z_init = max_z;
    min_z_init = min_z;
#endif


    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;
#ifdef THREEDIM
    T Lz = max_z - min_z;
#endif

    // safety_factor = 2 means we have a grid which has a grid point 2*dx from the boundary particle
    // Assuming a local approach, a grid point 2dx away from a particle will not influence this particle
    unsigned int safety_factor_x = sfx; // std::max(sfx, nonlocal_support);
    unsigned int safety_factor_y = sfy; // std::max(sfy, nonlocal_support);
    unsigned int safety_factor_z = sfz; // std::max(sfz, nonlocal_support);

    Nx = std::ceil(Lx * one_over_dx) + 1 + 2*safety_factor_x;
    Ny = std::ceil(Ly * one_over_dx) + 1 + 2*safety_factor_y;
#ifdef THREEDIM
    Nz = std::ceil(Lz * one_over_dx) + 1 + 2*safety_factor_z;
    grid_nodes = Nx*Ny*Nz;
#else
    grid_nodes = Nx*Ny;
#endif

    // save for remeshFixedCont
    Nx_init = Nx;
    Ny_init = Ny;
#ifdef THREEDIM
    Nz_init = Nz;
#endif

    low_x_init    = min_x - dx * safety_factor_x;
    high_x_init   = max_x + dx * safety_factor_x;
    low_y_init    = min_y - dx * safety_factor_y;
    high_y_init   = max_y + dx * safety_factor_y;
#ifdef THREEDIM
    T low_z_init  = min_z - dx * safety_factor_z;
    T high_z_init = max_z + dx * safety_factor_z;
#endif

#ifdef WARNINGS
    #ifdef THREEDIM
        debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );
    #else
        debug("               grid   = (", Nx, ", ", Ny, ")"  );
    #endif
#endif

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x_init, high_x_init, Nx);
    grid.y = linspace(low_y_init, high_y_init, Ny);
#ifdef THREEDIM
    grid.z = linspace(low_z_init, high_z_init, Nz);
#endif

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
#ifdef THREEDIM
    grid.zc = grid.z[0];
#endif

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}

void Simulation::remeshFixedCont(){

    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = (*min_x_it)(0);

    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);
#ifdef THREEDIM
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);
#endif

    T high_x;
    if (max_x < max_x_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_x_init-max_x)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +x reduction = ", reduction_factor);
#endif
        Nx     = Nx_init     - reduction_factor;
        high_x = high_x_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_x-max_x_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +x expansion = ", expansion_factor);
#endif
        Nx     = Nx_init     + expansion_factor;
        high_x = high_x_init + expansion_factor * dx;
    }

    T low_x;
    if (min_x > min_x_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_x-min_x_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -x reduction = ", reduction_factor);
#endif
        Nx    = Nx         - reduction_factor;
        low_x = low_x_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_x_init-min_x)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -x expansion = ", expansion_factor);
#endif
        Nx    = Nx         + expansion_factor;
        low_x = low_x_init - expansion_factor * dx;
    }

    T high_y;
    if (max_y < max_y_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_y_init-max_y)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +y reduction = ", reduction_factor);
#endif
        Ny     = Ny_init     - reduction_factor;
        high_y = high_y_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_y-max_y_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +y expansion = ", expansion_factor);
#endif
        Ny     = Ny_init     + expansion_factor;
        high_y = high_y_init + expansion_factor * dx;
    }

    T low_y;
    if (min_y > min_y_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_y-min_y_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -y reduction = ", reduction_factor);
#endif
        Ny    = Ny         - reduction_factor;
        low_y = low_y_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_y_init-min_y)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -y expansion = ", expansion_factor);
#endif
        Ny    = Ny         + expansion_factor;
        low_y = low_y_init - expansion_factor * dx;
    }

#ifdef THREEDIM
    T high_z;
    if (max_z < max_z_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_z_init-max_z)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +z reduction = ", reduction_factor);
#endif
        Nz     = Nz_init     - reduction_factor;
        high_z = high_z_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_z-max_z_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid +z expansion = ", expansion_factor);
#endif
        Nz     = Nz_init     + expansion_factor;
        high_z = high_z_init + expansion_factor * dx;
    }

    T low_z;
    if (min_z > min_z_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_z-min_z_init)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -z reduction = ", reduction_factor);
#endif
        Nz    = Nz         - reduction_factor;
        low_z = low_z_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_z_init-min_z)/dx - 1e-8*dx) );
#ifdef WARNINGS
        debug("               grid -z expansion = ", expansion_factor);
#endif
        Nz    = Nz         + expansion_factor;
        low_z = low_z_init - expansion_factor * dx;
    }
#endif



#ifdef THREEDIM
    grid_nodes = Nx*Ny*Nz;
#else
    grid_nodes = Nx*Ny;
#endif

#ifdef WARNINGS
    #ifdef THREEDIM
        debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );
    #else
        debug("               grid   = (", Nx, ", ", Ny, ")"  );
    #endif
#endif

    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x, high_x, Nx);
    grid.y = linspace(low_y, high_y, Ny);

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
#ifdef THREEDIM
    grid.zc = grid.z[0];
#endif

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

}
