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
    frame_dt = 1.0 / fps;

    lambda = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    mu     = E / (2.0*(1.0+nu));
    K = lambda + 2.0 * mu / dim;
    rho = density;

    wave_speed = std::sqrt(E/rho);
    dt_max = dt_max_coeff * dx / wave_speed;

    particles = Particles(Np);

#ifdef THREEDIM
    particle_volume = Lx*Ly*Lz / Np; // INITIAL particle volume V^0
#else
    particle_volume = Lx*Ly / Np;
#endif
    particle_mass = rho * particle_volume;

    nonlocal_l_sq = nonlocal_l * nonlocal_l;
    nonlocal_support = std::ceil(nonlocal_l / dx);

    mu_sqrt6 = mu * std::sqrt((T)6);

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

}

void Simulation::simulate(){

    std::cout << "------------------------------------------------ " << std::endl;
    std::cout << "---------------- This is Larsie ---------------- " << std::endl;
    std::cout << "------------------------------------------------ " << std::endl;

    createDirectory();

    // Write parameters to file for future reference
    std::ofstream infoFile(directory + sim_name + "/info.txt");
    infoFile << end_frame           << "\n"   // 0
             << fps                 << "\n"   // 1
             << dx                  << "\n"   // 2
             << mu                  << "\n"   // 3
             << lambda              << "\n"   // 4
             << vmin_factor         << "\n"   // 5
             << load_factor         << "\n";  // 6
    infoFile.close();

    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Precomputations
    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;

    // Precomputations for Drucker Prager
    // T sin_phi = std::sin(friction_angle / 180.0 * M_PI);
    // T alpha = std::sqrt(2.0/3.0) * 2.0 * sin_phi / (3.0 - sin_phi);
    // alpha_K_d_over_2mu = alpha * K * dim / (2*mu); // = alpha * bulk_modulus * dimension / (2*mu)
    // particles.cohesion_proj.resize(Np);  std::fill( particles.cohesion_proj.begin(),  particles.cohesion_proj.end(),  cohesion );

    // Lagrangian coordinates. Using assignment operator to copy
    particles.x0 = particles.x;

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;
    saveParticleData();
    while (frame < end_frame){
        std::cout << "Frame: " << frame << std::endl;
        std::cout << "               Step: " << current_time_step  << std::endl;
        std::cout << "               Time: " << time               << std::endl;
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

    // remesh();

    if (current_time_step == 0) {
        remeshFixedInit();
    } else {
        remeshFixedCont();
    }

    moveObjects();
    P2G();
    // calculateMassConservation();
    explicitEulerUpdate();
    // addExternalParticleGravity();
    G2P();
    deformationUpdate();
    // plasticity_projection();
    positionUpdate();
}



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
    explicitEulerUpdate_Optimized_Parallel(); // CURRENTLY NOT WORKING!
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

} // end updateDt





void Simulation::moveObjects(){
    for (InfinitePlate &obj : objects) {
        obj.move(dt, frame_dt, time);
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
    std::ofstream outFile(directory + sim_name + "/out_part_frame_" + extra + std::to_string(frame) + ".csv");

    outFile << "x"           << ","   // 0
            << "y"           << ","   // 1
            << "z"           << ","   // 2
            << "vx"          << ","   // 3
            << "vy"          << ","   // 4
            << "vz"          << ","   // 5
            << "pressure"    << ","   // 6
            << "devstress"   << ","   // 7
            << "eps_pl_vol"  << ","   // 8
            << "eps_pl_dev"  << ","   // 9
            << "eps_pl_dev_nonloc"  << ","   // 10
            << "delta_gamma"        << ","   // 11
            << "delta_gamma_nonloc" << "\n";   // 12
            // << "tau_xx"      << ","   // 13
            // << "tau_xy"      << ","   // 14
            // << "tau_yx"      << ","   // 15
            // << "tau_yy"      << ","   // 16
            // << "Fe_xx"       << ","   // 17
            // << "Fe_xy"       << ","   // 18
            // << "Fe_yx"       << ","   // 19
            // << "Fe_yy"       << "\n";   // 20

    TM I = TM::Identity();
    TM volavg_tau = TM::Zero();
    T Jsum = 0;
    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        TM tau; // particles.tau[p];
        if (elastic_model == NeoHookean)
            tau = NeoHookeanPiola(Fe) * Fe.transpose();
        else if (elastic_model == StvkWithHencky)
            tau = StvkWithHenckyPiola(Fe) * Fe.transpose();

        T J = Fe.determinant() * std::exp( particles.eps_pl_vol[p] );
        volavg_tau += tau * J;
        Jsum += J;

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * I;
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        outFile << particles.x[p](0)          << ","   // 0
                << particles.x[p](1)          << ","   // 1
            #ifdef THREEDIM
                << particles.x[p](2)          << ","   // 2
            #else
                << 0                          << ","
            #endif
                << particles.v[p](0)          << ","   // 3
                << particles.v[p](1)          << ","   // 4
            #ifdef THREEDIM
                << particles.v[p](2)          << ","   // 5
            #else
                << 0                          << ","
            #endif
                << pressure                   << ","   // 6
                << devstress                  << ","   // 7
                << particles.eps_pl_vol[p]    << ","   // 8
                << particles.eps_pl_dev[p]    << ","   // 9
                << particles.eps_pl_dev_nonloc[p]   << ","     // 10
                << particles.delta_gamma[p]         << ","     // 11
                << particles.delta_gamma_nonloc[p]  << "\n";   // 12
                // << tau(0,0)                   << ","   // 13
                // << tau(0,1)                   << ","   // 14
                // << tau(1,0)                   << ","   // 15
                // << tau(1,1)                   << ","   // 16
                // << Fe(0,0)                    << ","    // 17
                // << Fe(0,1)                    << ","    // 18
                // << Fe(1,0)                    << ","    // 19
                // << Fe(1,1)                    << "\n";  // 20
    } // end loop over particles
    outFile.close();

    volavg_tau /= Jsum;
    T volavg_pressure = -volavg_tau.trace() / dim;
    TM volavg_tau_dev = volavg_tau + volavg_pressure * I;
    T volavg_devstress = std::sqrt(3.0/2.0 * selfDoubleDot(volavg_tau_dev));
    std::ofstream outFile2(directory + sim_name + "/out_pq_frame_" + extra + std::to_string(frame) + ".csv");
    outFile2 << volavg_pressure    << ","
             << volavg_devstress   << "\n";
    outFile2.close();

    std::ofstream outFile3(directory + sim_name + "/last_written.txt");
    outFile3 << std::to_string(frame) << "\n";
    outFile3.close();

} // end saveParticleData()

void Simulation::saveGridData(std::string extra){
    std::ofstream outFile(directory + sim_name + "/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    outFile         << "x"       << ","
                    << "y"       << ","
                    << "z"       << ","
                    << "vx"      << ","
                    << "vy"      << ","
                    << "vz"      << ","
                    << "mass"    << "\n";
                //    << "delta_gamma" << ","; // OBS: Only if nonlocal strategy

#ifdef THREEDIM
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; k++){
                unsigned int index = ind(i,j,k);
                outFile << grid.x[i]               << "," // 0
                        << grid.y[j]               << "," // 1
                        << grid.z[k]               << "," // 2
                        << grid.v[index](0)        << "," // 3
                        << grid.v[index](1)        << "," // 4
                        << grid.v[index](2)        << "," // 5
                        << grid.mass[index]        << "\n";   // 6
                    //    << grid.delta_gamma[index] << ",";
            } // end for k
        } // end for j
    } // end for i
    outFile.close();
#else
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            unsigned int index = ind(i,j);
            outFile << grid.x[i]               << "," // 0
                    << grid.y[j]               << "," // 1
                    << 0                       << ","
                    << grid.v[index](0)        << "," // 3
                    << grid.v[index](1)        << "," // 4
                    << 0                       << ","
                    << grid.mass[index]        << "\n";  // 6
                //    << grid.delta_gamma[index] << ","; // 7
        } // end for j
    } // end for i
    outFile.close();
#endif

} // end saveGridData()


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
