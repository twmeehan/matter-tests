// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

Simulation::Simulation(){
    current_time_step = 0;
    time  = 0;
    frame = 0;
    exit  = 0;

    is_initialized = false;

    runtime_p2g = 0;
    runtime_g2p = 0;
    runtime_euler = 0;
    runtime_defgrad = 0;

    // create grid
    grid = Grid();
}


void Simulation::initialize(bool save, std::string dir, std::string name){

    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "    88b           d88                                                              " << std::endl;
    std::cout << "    888b         d888                  aa          aa                              " << std::endl;
    std::cout << "    88`8b       d8'88                  88          88                              " << std::endl;
    std::cout << "    88 `8b     d8' 88  ,adPPYYba,  aaaa88aaaa  aaaa88aaaa   ,adPPYba,  8b,dPPYba,  " << std::endl;
    std::cout << "    88  `8b   d8'  88  aa     `Y8  aaaa88aaaa  8888888888  a8P_____88  88P     Y8  " << std::endl;
    std::cout << "    88   `8b d8'   88  ,adPPPPP88      88          88      adPPPPP88   88          " << std::endl;
    std::cout << "    88    `888'    88  88,    ,88      aa          aa      a8b         88          " << std::endl;
    std::cout << "    88     `8'     88   `adPPYba,                           `adPPYba   88          " << std::endl;
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;

    save_sim = save;
    directory = dir;
    sim_name = name;

    if (save_sim)
        createDirectory();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;

    is_initialized = true;
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

// NB: Simulation::initialize(...) must be called before Simulation::simulate()
void Simulation::simulate(){

    if (!is_initialized){
        debug("Simulation not initialized. Call the initialize(...) function in your input file.");
        return;
    }

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

    lambda = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) ); // first Lame parameter
    mu = E / (2.0*(1.0+nu)); // shear modulus
    K = calculateBulkModulus(); // bulk modulus
    wave_speed = std::sqrt(E/rho); // elastic wave speed

    dt_max = cfl_elastic * dx / wave_speed;

    frame_dt = 1.0 / fps;

    gravity_final = gravity;

    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;

    if (use_von_mises_q){
        q_prefac = sqrt3/sqrt2;
        d_prefac = sqrt2/sqrt3;
    } else {
        q_prefac = 1.0/sqrt2;
        d_prefac = sqrt2;
    }
    e_mu_prefac = 2*q_prefac          * mu;
    f_mu_prefac = 2*q_prefac/d_prefac * mu;
    rma_prefac  = 2*q_prefac*q_prefac;

    fac_Q = I_ref / (grain_diameter*std::sqrt(rho_s));

    if (use_material_friction)
        std::fill(particles.muI.begin(), particles.muI.end(), mu_1);

    debug("Number of particles: ", Np);
    debug("Grid spacing dx:     ", dx);
    debug("Elastic wave speed:  ", wave_speed);
    debug("Maximum dt:          ", dt_max);
    debug("Particle volume:     ", particle_volume);
    debug("Particle mass:       ", particle_mass);

    // particles.x0 = particles.x;  // Lagrangian coordinates, using assignment operator to copy

    // nonlocal_l_sq = nonlocal_l * nonlocal_l;
    // nonlocal_support = std::ceil(nonlocal_l / dx);

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;

    if (save_sim){
        saveInfo();
        saveParticleData();
    }

    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    while (frame < end_frame){
        if (!reduce_verbose){
            std::cout << "Frame: "               << frame  << " / "    << end_frame  << std::endl;
            std::cout << "               Name: " << sim_name           << std::endl;
            std::cout << "               Step: " << current_time_step  << std::endl;
            std::cout << "               Time: " << time   << " -> "   << (frame+1)*frame_dt << std::endl;
        }
        advanceStep();
        if (exit == 1)
            return;
        time += dt;
        current_time_step++;
        if( time - frame_dt*(frame+1) >= 0 ){
            frame++;
            std::cout << "End of frame " << frame << std::endl;
            if (save_sim){
                saveParticleData();
                if (save_grid)
                    saveGridData();
                // saveAvgData();
            }

        }
        if (std::abs(final_time-time) < 1e-15 || final_time < time){
            std::cout << "The simulation ended at time = " << time << std::endl;
            if (save_sim){
                saveParticleData();
                if (save_grid)
                    saveGridData();
                // saveAvgData();
            }
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
            if (delete_last_particle)
                deleteLastParticle(delete_last_particle);
    }else{ //  not pbc
        if (current_time_step == 0) {
            remeshFixedInit(2,2,2);
            if (delete_last_particle)
                deleteLastParticle(delete_last_particle);
        } else {
            remeshFixedCont();
        }
    }

    moveObjects();

    if (pbc){
        PBCAddParticles(4);
    }

    timer t_p2g; t_p2g.start();
    P2G();
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();


    // checkMassConservation();
    // checkMomentumConservation();

    timer t_euler; t_euler.start();
    explicitEulerUpdate();
    t_euler.stop(); runtime_euler += t_euler.get_timing();

    // addExternalParticleGravity();

    if (pbc){
        PBCDelParticles();
    }


    timer t_g2p; t_g2p.start();
    G2P();
    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();

    if (musl == true)
        MUSL();

    timer t_defgrad; t_defgrad.start();
    deformationUpdate();
    t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();

    // plasticity_projection(); // if nonlocal approach (deprecated)

    positionUpdate();

} // end advanceStep



void Simulation::moveObjects(){
    for (ObjectPlate &obj : plates) {
        obj.move(dt, frame_dt, time);
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA HELPER FUNCTIONS //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

T Simulation::calculateBulkModulus(){
    return lambda + 2.0 * mu / dim;
}

void Simulation::checkMomentumConservation(){
    TV momentum_grid = TV::Zero();
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
        #ifdef THREEDIM
            for(int k=0; k<Nz; k++)
                momentum_grid += grid.mass[ind(i,j,k)] * grid.v[ind(i,j,k)];
        #else
            momentum_grid += grid.mass[ind(i,j)] * grid.v[ind(i,j)];
        #endif
        }
    }

    TV momentum_particle = TV::Zero();
    for(int p=0; p<Np; p++){
        momentum_particle += particle_mass * particles.v[p];
    }
    debug("               Total part momentum = ", momentum_grid.norm());
    debug("               Total grid momentum = ", momentum_particle.norm());

    if ( (momentum_grid-momentum_particle).norm() > 1e-5 * momentum_particle.norm() ){
        debug("MOMENTUM NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}

void Simulation::checkMassConservation(){
    T particle_mass_total = particle_mass*Np;
    T grid_mass_total = 0;
    for(auto&& m: grid.mass)
        grid_mass_total += m;

    debug("               Total grid mass = ", grid_mass_total    );
    debug("               Total part mass = ", particle_mass_total);

    if ( std::abs(grid_mass_total-particle_mass_total) > 1e-5 * particle_mass_total ){
        debug("MASS NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}


void Simulation::deleteLastParticle(unsigned int n){

    for (unsigned int i = 0; i < n; i++) {

        auto new_part_x = particles.x;
        new_part_x.pop_back();

        Np -= 1;

        if (new_part_x.size() != Np){
            debug("PARTICLE NUMBER MISMATCH!!!");
            exit = 1;
            return;
        }

        particles = Particles(Np);
        particles.x = new_part_x;

    }
}


// This function is to be used in explicitEulerUpdate after boundaryCollision
// It must be hard-coded to choice
void Simulation::overwriteGridVelocity(TV Xi, TV& vi){
    T y_start = Ly - 0.25*dx;
    T width = 2*dx;
    T v_imp = 0.1; // positive value means tension
    if (Xi(1) > y_start - width + v_imp * time)
        vi(1) = v_imp;
    if (Xi(1) < 0       + width - v_imp * time)
        vi(1) = -v_imp;
}


// This function must be hard-coded to choice
void Simulation::addExternalParticleGravity(){
    // // 1. Transfer grid velocity to particles
    // G2P();
    // // 2. Apply gravity on particle velocity
    // for(int p=0; p<Np; p++)
    //     particles.v[p] += dt * (-2.0*amplitude*particles.x0[p]);
    // // 3. Transfer particle velocity back to grid
    // P2G();
} // end addExternalParticleGravity


// // TODO:
// void Simulation::boundaryCorrection(T xi, T yi, T& vxi, T& vyi){
//
//     // trial step
//     T x_next = xi + vxi * dt;
//     T y_next = yi + vyi * dt;
//     moveObjects();
//
//     for (ObjectPlate &obj : objects_plate) {
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
