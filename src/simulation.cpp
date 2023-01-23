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

    // create grid
    grid = Grid();
}


void Simulation::initialize(T E, T nu, T density){

    lambda = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    mu     = E / (2.0*(1.0+nu));
    K = lambda + 2.0 * mu / dim;
    rho = density;
    wave_speed = std::sqrt(E/rho);

}


void Simulation::createDirectory(){

    if (mkdir(directory.c_str(), 0777) == -1)
        std::cerr << "Directory " << directory << " has already been created:  " << strerror(errno) << std::endl;
    else
        std::cout << "Directory " << directory << " was created now"<< std::endl;
    if (mkdir((directory + sim_name).c_str(), 0777) == -1)
        std::cerr << "Simulation " << sim_name << " has already been created: " << strerror(errno) << std::endl;
    else
        std::cout << "Simulation " << sim_name << " was created now" << std::endl;

    std::string in  = "../src/mpm.cpp";
    std::string out = (directory + sim_name + "/initial_data.cpp");
    bool check = copy_file(in, out);
    if (!check){
        std::cerr << "Initial data " << in << " was NOT successfully copied to " << out << std::endl;
    }

}

void Simulation::simulate(){

    #if DIMENSION == 3
        debug("This is a 3D simulation.");
    #elif DIMENSION == 2
        debug("This is a 2D simulation.");
    #else
        #error Unsupported spline degree
    #endif

    #if SPLINEDEG == 3
      apicDinverse = 3.0/(dx*dx);
      debug("Using cubic splines.");
    #elif SPLINEDEG == 2
      apicDinverse = 4.0/(dx*dx);
      debug("Using quadratic splines.");
    #elif SPLINEDEG == 1
        apicDinverse = 0; // NB not implemented
        debug("Using linear hat functions.");
    #else
        #error Unsupported spline degree
    #endif

    // Precomputations
    frame_dt = 1.0 / fps;

    gravity_final = gravity;

    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;

    nonlocal_l_sq = nonlocal_l * nonlocal_l;
    nonlocal_support = std::ceil(nonlocal_l / dx);

    fac_Q = in_numb_ref / (grain_diameter*std::sqrt(rho_s)); // NB: Use 2 * grain diameter if using the other definiton

    debug("Num of particles = ", Np);
    debug("dx               = ", dx);
    debug("Wave speed       = ", wave_speed);
    debug("dt_max           = ", dt_max);
    debug("particle_volume  = ", particle_volume);
    debug("particle_mass    = ", particle_mass);

    std::cout << "------------------------------------------------ " << std::endl;
    std::cout << "---------------- This is Larsie ---------------- " << std::endl;
    std::cout << "------------------------------------------------ " << std::endl;

    createDirectory();

    saveInfo();

    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Lagrangian coordinates. Using assignment operator to copy
    particles.x0 = particles.x;

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;
    saveParticleData();
    while (frame < end_frame){
        std::cout << "Frame: "               << frame  << " / "    << end_frame  << std::endl;
        std::cout << "               Name: " << sim_name           << std::endl;
        std::cout << "               Step: " << current_time_step  << std::endl;
        std::cout << "               Time: " << time   << " -> "   << (frame+1)*frame_dt << std::endl;
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
            saveParticleData();
            saveGridData();
            break;
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Simulation took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " milliseconds" << std::endl;
    debug("Runtime P2G     = ", runtime_p2g     * 1000.0, " milliseconds");
    debug("Runtime G2P     = ", runtime_g2p     * 1000.0, " milliseconds");
    debug("Runtime Euler   = ", runtime_euler   * 1000.0, " milliseconds");
    debug("Runtime DefGrad = ", runtime_defgrad * 1000.0, " milliseconds");
}





void Simulation::advanceStep(){
    updateDt();

    if (pbc){
        if (current_time_step == 0)
            remeshFixed(4);
    }else{
        if (current_time_step == 0) {
            remeshFixedInit(2,2,2);
        } else {
            remeshFixedCont();
        }
    }

    moveObjects();

    if (pbc){
        // PBCAddParticles1D();
        PBCAddParticles(4);
    }

    P2G();
    // calculateMassConservation();
    explicitEulerUpdate();
    // addExternalParticleGravity();

    if (pbc){
        PBCDelParticles();
    }

    G2P();
    deformationUpdate();
    // plasticity_projection();   // if nonlocal approach

    positionUpdate();

} // end advanceStep



void Simulation::P2G(){
    timer t_p2g; t_p2g.start();
    // P2G_Baseline();
    // P2G_Optimized();
    P2G_Optimized_Parallel();
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();
} // end P2G

void Simulation::explicitEulerUpdate(){
    timer t_euler; t_euler.start();
    // explicitEulerUpdate_Baseline();
    // explicitEulerUpdate_Optimized();
    explicitEulerUpdate_Optimized_Parallel();
    t_euler.stop(); runtime_euler += t_euler.get_timing();
}

void Simulation::deformationUpdate(){
    timer t_defgrad; t_defgrad.start();
    // deformationUpdate_Baseline();
    deformationUpdate_Parallel();
    t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();
}

void Simulation::G2P(){
    timer t_g2p; t_g2p.start();
    // G2P_Baseline();
    // G2P_Optimized();
    G2P_Optimized_Parallel();
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

    if (max_speed >= wave_speed){
        debug("DETECTED SPEED LARGER THAN ELASTIC WAVE SPEED!!!");
        exit = 1;
        return;
    }

#ifdef WARNINGS
    debug("               dt_max = ", dt_max);
#endif

    if (std::abs(max_speed) > 1e-10){
        T dt_cfl = cfl * dx / max_speed;
#ifdef WARNINGS
        debug("               dt_cfl = ", dt_cfl);
#endif
        dt = std::min(dt_cfl, dt_max);
    } else {
        dt = dt_max;
#ifdef WARNINGS
        debug("               dt_cfl = not computed, max_speed too low");
#endif
    }

    dt = std::min(dt, frame_dt*(frame+1) - time);
    dt = std::min(dt, final_time         - time);

#ifdef WARNINGS
    debug("               dt     = ", dt    );
#endif

    // if (dt > dt_cfl){
    //     debug("TIME STEP IS TOO BIG COMPARED TO CFL!!!");
    //     exit = 1;
    //     return;
    // }
    // if (dt > dt_max){
    //     debug("TIME STEP IS TOO BIG COMPARED TO ELASTIC WAVE SPEED!!!");
    //     exit = 1;
    //     return;
    // }

    if (time < gravity_time){
        gravity = gravity_final * time/gravity_time;
    }
    else{
        gravity = gravity_final;
    }

} // end updateDt



void Simulation::moveObjects(){
    for (InfinitePlate &obj : objects) {
        obj.move(dt, frame_dt, time);
    }
}

void Simulation::positionUpdate(){

    #pragma omp parallel for num_threads(n_threads)
    for(int p=0; p<Np; p++){

        //// Position is updated according to PIC velocities
        particles.x[p] = particles.x[p] + dt * particles.pic[p];

        //// Velicity is updated
        if (flip_ratio < -1){ // APIC
            particles.v[p] = particles.pic[p];
        } else if (flip_ratio < 0){ // AFLIP
            particles.v[p] = (-flip_ratio) * ( particles.v[p] + particles.flip[p] ) + (1 - (-flip_ratio)) * particles.pic[p];
        } else{ // PIC-FLIP
            particles.v[p] =   flip_ratio  * ( particles.v[p] + particles.flip[p] ) + (1 -   flip_ratio)  * particles.pic[p];
        }


        if (pbc){
            if (particles.x[p](0) > Lx){
                particles.x[p](0) = particles.x[p](0) - Lx;
            }
            else if (particles.x[p](0) < 0){
                particles.x[p](0) = Lx + particles.x[p](0);
            }
        }

        if (pbc_special){
            // if (particles.x[p](0) > 0.8){
            //     particles.x[p](0)  = -0.1;
            //     particles.x[p](1) +=  0.27;
            // }
            if (particles.x[p](0) > 1.8){
                particles.x[p](0)  = -0.18;
                particles.x[p](1) +=  0.5;

                particles.v[p](0) *= 0.1;
            }
        }

    } // end loop over particles
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA FUNCTIONS //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// // Currently not working!!
// void Simulation::boundaryCorrection(T xi, T yi, T& vxi, T& vyi){
//
//     // trial step
//     T x_next = xi + vxi * dt;
//     T y_next = yi + vyi * dt;
//     moveObjects();
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
#ifdef THREEDIM
    void Simulation::overwriteGridVelocity(T xi, T yi, T zi, TV& vi){
#else
    void Simulation::overwriteGridVelocity(T xi, T yi, TV& vi){
#endif
    T y_start = Ly - 0.25*dx;
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
        #ifdef THREEDIM
            for(int k=0; k<Nz; k++)
                momentum += grid.mass[ind(i,j,k)] * grid.v[ind(i,j,k)];
        #else
            momentum += grid.mass[ind(i,j)] * grid.v[ind(i,j)];
        #endif

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
