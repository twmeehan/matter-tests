#include "simulation.hpp"

Simulation::Simulation(){
    current_time_step = 0;
    time  = 0;
    frame = 0;
    exit  = 0;

    runtime_p2g = 0;
    runtime_g2p = 0;
    runtime_euler = 0;
    runtime_defgrad = 0;

    // Default values
    Nx = 3;
    Ny = 3;
    Nz = 3;
    grid = Grid(Nx, Ny, Nz);
}


void Simulation::initialize(T E, T nu, T density){

    dim = 3;

    lambda = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    mu     = E / (2.0*(1.0+nu));
    K = lambda + 2.0 * mu / dim;
    rho = density;

    wave_speed = std::sqrt(E/rho);
    dt_max = 0.1 * dx / wave_speed;

    particle_volume = std::pow(L, dim) / Np; // INITIAL particle volume V^0
    particle_mass = rho * particle_volume;

    particles = Particles(Np);


}




void Simulation::simulate(){

    std::cout << "------------------------------------------------ " << std::endl;
    std::cout << "---------------- This is Larsie ---------------- " << std::endl;
    std::cout << "------------------------------------------------ " << std::endl;

    // Create Directories:
    if (mkdir("dumps", 0777) == -1)
        std::cerr << "Directory dumps already created:  " << strerror(errno) << std::endl;
    else
        std::cout << "directory dumps created" << std::endl;
    if (mkdir(("dumps/" + sim_name).c_str(), 0777) == -1)
        std::cerr << "Directory " << sim_name << " already created: " << strerror(errno) << std::endl;
    else
        std::cout << "Directory " << sim_name << " created" << std::endl;

    // Write parameters to file for future reference
    int emodel, pmodel;
    if (elastic_model == StvkWithHencky)
        emodel = 1;
    else if (elastic_model == NeoHookean)
        emodel = 2;
    else
        emodel = -1;
    if (plastic_model == VonMises)
        pmodel = 1;
    else if (plastic_model == NoPlasticity)
        pmodel = 2;
    else
        pmodel = -1;
    std::ofstream infoFile("dumps/" + sim_name + "/info.txt");
    infoFile << end_frame           << "\n"  // 0
             << frame_dt            << "\n"  // 1
             << dx                  << "\n"  // 2
             << mu                  << "\n"  // 3
             << lambda              << "\n"  // 4
             << emodel              << "\n"  // 5
             << pmodel              << "\n"  // 6
             << xi                  << "\n"  // 7
             << yield_stress_orig   << "\n"  // 8
             << yield_stress_min    << "\n"  // 9
             << friction_angle      << "\n"  // 10
             << cohesion            << "\n"  // 11
             << reg_length          << "\n"  // 12
             << reg_const           << "\n"; // 13
    infoFile.close();

    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Precomputations
    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;
    reg_const_length_sq = reg_const * reg_length * reg_length;

    // Precomputations for Drucker Prager
    T sin_phi = std::sin(friction_angle / 180.0 * M_PI);
    T alpha = std::sqrt(2.0/3.0) * 2.0 * sin_phi / (3.0 - sin_phi);
    alpha_K_d_over_2mu = alpha * K * dim / (2*mu); // = alpha * bulk_modulus * dimension / (2*mu)
    particles.cohesion_proj.resize(Np);  std::fill( particles.cohesion_proj.begin(),  particles.cohesion_proj.end(),  cohesion );

    // Lagrangian coordinates. Using assignment operator to copy
    particles.x0 = particles.x;

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;
    saveParticleData();
    saveGridData();
    while (frame < end_frame){
        std::cout << "Step: " << current_time_step << "    Time: " << time << std::endl;
        advanceStep();
        if (exit == 1)
            return;
        time += dt;
        current_time_step++;
        if( std::abs(time - frame_dt*(frame+1)) < 1e-15 ){
            frame++;
            std::cout << "Saving frame " << frame << std::endl;
            saveParticleData();
            saveGridData();
        }
        if (std::abs(final_time-time) < 1e-15 || final_time < time){
            std::cout << "The simulation ended successfully at time = " << time << std::endl;
            break;
        }
    }
    saveParticleData();
    saveGridData();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Simulation took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " milliseconds" << std::endl;
    debug("Runtime P2G     = ", runtime_p2g     * 1000.0, " milliseconds");
    debug("Runtime G2P     = ", runtime_g2p     * 1000.0, " milliseconds");
    debug("Runtime Euler   = ", runtime_euler   * 1000.0, " milliseconds");
    debug("Runtime DefGrad = ", runtime_defgrad * 1000.0, " milliseconds");
}





void Simulation::advanceStep(){
    updateDt();
    remesh();
    moveObjects(dt);
    P2G();
    calculateMassConservation();
    explicitEulerUpdate();
    // addExternalParticleGravity();
    G2P();              // due to "regularization", G2P must come before deformationUpdate!!!
    deformationUpdate();
    positionUpdate();
}



void Simulation::P2G(){
    timer t_p2g; t_p2g.start();
    // P2G_Baseline();
    P2G_Optimized();
    // P2G_Optimized_Parallel();
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();
} // end P2G

void Simulation::explicitEulerUpdate(){
    timer t_euler; t_euler.start();
    // explicitEulerUpdate_Baseline();
    explicitEulerUpdate_Optimized();
    // explicitEulerUpdate_Optimized_Parallel(); // CURRENTLY NOT WORKING!
    t_euler.stop(); runtime_euler += t_euler.get_timing();
}

void Simulation::deformationUpdate(){
    timer t_defgrad; t_defgrad.start();
    deformationUpdate_Baseline();
    // deformationUpdate_Parallel();
    t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();
}

void Simulation::G2P(){
    timer t_g2p; t_g2p.start();
    // G2P_Baseline();
    G2P_Optimized();
    // G2P_Optimized_Parallel();
    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();
}






void Simulation::updateDt(){

    // T max_speed = std::sqrt((particles.vx.array().square() + particles.vy.array().square()).maxCoeff());
    auto max_velocity_it = std::max_element( particles.v.begin(), particles.v.end(),
                                             []( const TV &v1, const TV &v2 )
                                             {
                                                 return v1.squaredNorm() < v2.squaredNorm();
                                             } );
    T max_speed = (*max_velocity_it).norm();
    // ///////////// DEBUG /////////////
    // debug("max_speed (iterator) = ", max_speed);
    // max_speed = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.v[p].norm() > max_speed){
    //         max_speed = particles.v[p].norm();
    //     }
    // }
    // debug("max_speed (debug)    = ", max_speed);
    // /////////////////////////////////

    T dt_cfl = cfl * dx / max_speed;
    dt = std::min(dt_cfl, dt_max);
    dt = std::min(dt, frame_dt*(frame+1) - time);
    dt = std::min(dt, final_time         - time);
    debug("               dt_cfl = ", dt_cfl);
    debug("               dt_max = ", dt_max);
    debug("               dt     = ", dt    );

    if (dt > dt_cfl){
        debug("TIME STEP IS TOO BIG COMPARED TO CFL!!!");
        exit = 1;
        return;
    }
    if (dt > dt_max){
        debug("TIME STEP IS TOO BIG COMPARED TO ELASTIC WAVE SPEED!!!");
        exit = 1;
        return;
    }
    if (max_speed >= wave_speed){
        debug("DETECTED SPEED LARGER THAN ELASTIC WAVE SPEED!!!");
        exit = 1;
        return;
    }
} // end updateDt





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

    // ///////////// DEBUG /////////////
    // debug("max_x (iterator) = ", max_x);
    // debug("max_y (iterator) = ", max_y);
    // debug("min_x (iterator) = ", min_x);
    // debug("min_y (iterator) = ", min_y);
    //
    // max_x = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](0) > max_x){
    //         max_x = particles.x[p](0);
    //     }
    // }
    // debug("max_x (debug)    = ", max_x);
    //
    // max_y = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](1) > max_y){
    //         max_y = particles.x[p](1);
    //     }
    // }
    // debug("max_y (debug)    = ", max_y);
    //
    // min_x = 1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](0) < min_x){
    //         min_x = particles.x[p](0);
    //     }
    // }
    // debug("min_x (debug)    = ", min_x);
    //
    // min_y = 1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](1) < min_y){
    //         min_y = particles.x[p](1);
    //     }
    // }
    // debug("min_y (debug)    = ", min_y);
    // /////////////////////////////////

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;
    T Lz = max_z - min_z;

    // ACTUAL midpoints of particles
    T mid_x = min_x + Lx / 2.0;
    T mid_y = min_y + Ly / 2.0;
    T mid_z = min_z + Lz / 2.0;

    // NEW number of grid-dx's per side length
    unsigned int safety_factor = 3;
    Nx = std::ceil(Lx * one_over_dx) + safety_factor;
    Ny = std::ceil(Ly * one_over_dx) + safety_factor;
    Nz = std::ceil(Lz * one_over_dx) + safety_factor;

    T low_x  = mid_x - Nx*dx/2.0;
    T low_y  = mid_y - Ny*dx/2.0;
    T low_z  = mid_z - Nz*dx/2.0;
    T high_x = mid_x + Nx*dx/2.0;
    T high_y = mid_y + Ny*dx/2.0;
    T high_z = mid_z + Nz*dx/2.0;

    // NEW grid dimensions N = (L / dx + 1)
    Nx++;
    Ny++;
    Nz++;

    debug(  "               grid   = (", Nx, ", ", Ny, ", ", Nz, ")"  );

    // Linspace x and y
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


void Simulation::moveObjects(T delta_t){
    for (InfinitePlate &obj : objects) {
        obj.move(delta_t);
    }
}



void Simulation::positionUpdate(){
    for(int p=0; p<Np; p++){

        // Position is updated according to PIC velocities
        particles.x[p] = particles.x[p] + dt * particles.pic[p];

        // New particle velocity is a FLIP-PIC combination
        particles.v[p] = flip_ratio * ( particles.v[p] + particles.flip[p] ) + (1 - flip_ratio) * particles.pic[p];

    } // end loop over particles
}





void Simulation::saveParticleData(std::string extra){
    std::ofstream outFile("dumps/" + sim_name + "/out_part_frame_" + extra + std::to_string(frame) + ".csv");

    outFile << "x"           << ","   // 0
            << "y"           << ","   // 1
            << "z"           << ","   // 2
            << "vx"          << ","   // 3
            << "vy"          << ","   // 4
            << "vz"          << ","   // 5
            << "eps_pl_dev"  << ","   // 6
            << "eps_pl_vol"  << ","   // 7
            << "reg"         << ","   // 8
            << "pressure"    << ","   // 9
            << "devstress"   << ","   // 10
            << "tau_xx"      << ","   // 11
            << "tau_xy"      << ","   // 12
            << "tau_xz"      << ","   // 13
            << "tau_yx"      << ","   // 14
            << "tau_yy"      << ","   // 15
            << "tau_yz"      << ","   // 16
            << "tau_zx"      << ","   // 17
            << "tau_zy"      << ","   // 18
            << "tau_zz"      << ","   // 19
            << "Fe_xx"       << ","   // 20
            << "Fe_xy"       << ","   // 21
            << "Fe_xz"       << ","   // 22
            << "Fe_yx"       << ","   // 23
            << "Fe_yy"       << ","   // 24
            << "Fe_yz"       << ","   // 25
            << "Fe_zx"       << ","   // 26
            << "Fe_zy"       << ","   // 27
            << "Fe_zz"       << "\n"; // 28

    for(int p = 0; p < Np; p++){

        TM I = TM::Identity();

        TM Fe = particles.F[p];

        TM tau; // particles.tau[p];
        if (elastic_model == NeoHookean)
            tau = NeoHookeanPiola(Fe) * Fe.transpose();
        else if (elastic_model == StvkWithHencky)
            tau = StvkWithHenckyPiola(Fe) * Fe.transpose();

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * I;
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        T reg = reg_length * reg_length * particles.regularization[p];

        outFile << particles.x[p](0)          << ","   // 0
                << particles.x[p](1)          << ","   // 1
                << particles.x[p](2)          << ","   // 2
                << particles.v[p](0)          << ","   // 3
                << particles.v[p](1)          << ","   // 4
                << particles.v[p](2)          << ","   // 5
                << particles.eps_pl_dev[p]    << ","   // 6
                << particles.eps_pl_vol[p]    << ","   // 7
                << reg                        << ","   // 8
                << pressure                   << ","   // 9
                << devstress                  << ","   // 10
                << tau(0,0)                   << ","   // 11
                << tau(0,1)                   << ","   // 12
                << tau(0,2)                   << ","   // 13
                << tau(1,0)                   << ","   // 14
                << tau(1,1)                   << ","   // 15
                << tau(1,2)                   << ","   // 16
                << tau(2,0)                   << ","   // 17
                << tau(2,1)                   << ","   // 18
                << tau(2,2)                   << ","   // 19
                << Fe(0,0)                    << ","    // 20
                << Fe(0,1)                    << ","    // 21
                << Fe(0,2)                    << ","    // 22
                << Fe(1,0)                    << ","    // 23
                << Fe(1,1)                    << ","    // 24
                << Fe(1,2)                    << ","    // 25
                << Fe(2,0)                    << ","    // 26
                << Fe(2,1)                    << ","    // 27
                << Fe(2,2)                    << "\n";  // 28
    }
}

void Simulation::saveGridData(std::string extra){
    std::ofstream outFile("dumps/" + sim_name + "/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    outFile         << "x"       << ","
                    << "y"       << ","
                    << "z"       << ","
                    << "vx"      << ","
                    << "vy"      << ","
                    << "vz"      << ","
                    << "mass"    << ","
                    << "reg"     << "\n";

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; k++){
                outFile << grid.x[i]             << ","
                        << grid.y[j]             << ","
                        << grid.z[k]             << ","
                        << grid.v[ind(i,j,k)](0) << ","
                        << grid.v[ind(i,j,k)](1) << ","
                        << grid.v[ind(i,j,k)](2) << ","
                        << grid.mass[ind(i,j,k)] << ","
                        << grid.regularization[ind(i,j,k)] << "\n";
            }
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA FUNCTIONS //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Simulation::boundaryCollision(T xi, T yi, T zi, TV& vi){

    // TODO: Make this a vector-based function

    // Reference to velocity components
    T& vxi = vi(0);
    T& vyi = vi(1);
    T& vzi = vi(2);

    // New positions
    xi += dt * vxi;
    yi += dt * vyi;
    zi += dt * vzi;

    for (InfinitePlate &obj : objects) {
        bool colliding = obj.inside(xi, yi);
        if (colliding) {
            T vx_rel = vxi - obj.vx_object;
            T vy_rel = vyi - obj.vy_object;

            if (boundary_condition == STICKY) {
                vx_rel = 0;
                vy_rel = 0;
            } // end STICKY

            else if (boundary_condition == SLIP) {
                if (obj.plate_type == top || obj.plate_type == bottom){
                    // tangential velocity is the x component
                    // normal component must be set to zero
                    vy_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vy_rel);
                        else
                            vx_rel += friction * std::abs(vy_rel);
                    }
                } // end top and bottom plate
                else if (obj.plate_type == left || obj.plate_type == right){
                    // tangential velocity is the y comp
                    // normal component must be set to zero
                    vx_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vx_rel);
                        else
                            vy_rel += friction * std::abs(vx_rel);
                    }
                } // end left or right plate
            } // end SLIP

            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vxi = vx_rel + obj.vx_object;
            vyi = vy_rel + obj.vy_object;
        } // end if colliding

    } // end iterator over objects

} // end boundaryCollision


// // Currently not working!!
// void Simulation::boundaryCorrection(T xi, T yi, T& vxi, T& vyi){
//
//     // trial step
//     T x_next = xi + vxi * dt;
//     T y_next = yi + vyi * dt;
//     moveObjects(dt);
//
//     for (InfinitePlate &obj : objects) {
//         bool colliding = obj.inside(x_next, y_next);
//         if (colliding) {
//             if (obj.plate_type == upper ){
//                 debug(obj.name, " dist = ", obj.distance(x_next, y_next));
//                 vyi += obj.distance(x_next, y_next) / dt; // distance is negative since grid point is inside object
//             }
//             else if (obj.plate_type == lower){
//                 debug(obj.name, " dist = ", obj.distance(x_next, y_next));
//                 vyi -= obj.distance(x_next, y_next) / dt;
//             }
//         } // end if colliding
//
//     } // end iterator over objects
//
//     // Correct back
//     moveObjects(-dt);
//
// } // end boundaryCorrection




// This function is to be used in explicitEulerUpdate after boundaryCollision
void Simulation::overwriteGridVelocity(T xi, T yi, T zi, TV& vi){
    T y_start = L - 0.25*dx;
    T width = 2*dx;
    T v_imp = 0.1; // positive value means tension
    if (yi > y_start - width + v_imp * time)
        vi(1) = v_imp;
    if (yi < 0       + width - v_imp * time)
        vi(1) = -v_imp;
}

void Simulation::addExternalParticleGravity(){
    // 1. Transfer grid velocity to particles
    G2P();
    // 2. Apply gravity on particle velocity
    for(int p=0; p<Np; p++)
        particles.v[p] += dt * (-2.0*amplitude*particles.x0[p]);
    // 3. Transfer particle velocity back to grid
    P2G();
} // end addExternalParticleGravity


// These functions are only for validating conservation laws
void Simulation::calculateMomentumOnParticles(){
    TV momentum = TV::Zero();
    for(int p=0; p<Np; p++){
        momentum += particle_mass * particles.v[p];
    }
    debug("               total part momentum = ", momentum.norm());
}
void Simulation::calculateMomentumOnGrid(){
    TV momentum = TV::Zero();
    for(int i=0; i<Nx; i++)
        for(int j=0; j<Ny; j++)
            for(int k=0; k<Nz; k++)
                momentum += grid.mass[ind(i,j,k)] * grid.v[ind(i,j,k)];
    debug("               total grid momentum = ", momentum.norm());
}

void Simulation::calculateMassConservation(){
    T particle_mass_total = particle_mass*Np;
    T grid_mass_total = 0;
    for(auto&& m: grid.mass)
        grid_mass_total += m;
    debug("               total grid mass = ", grid_mass_total    );
    debug("               total part mass = ", particle_mass_total);
    if ( std::abs(grid_mass_total-particle_mass_total) > 1e-5 * particle_mass_total ){
        debug("MASS NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}
