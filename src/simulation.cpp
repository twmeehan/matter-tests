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

    /////////// Default - will be overwritten by remesh if used ////////////
    Nx = 21; // N = (L / dx + 1)
    Ny = 21;
    grid = Grid(Nx, Ny);
    /////////////////////////////////////////////////////////////////////////
}


void Simulation::initialize(T E, T nu, T density){

    std::cout << "------------------------------------------------ " << std::endl;
    std::cout << "---------------- This is Larsie ---------------- " << std::endl;
    std::cout << "------------------------------------------------ " << std::endl;

    mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    lambda = E / (2.0*(1.0+nu));
    rho = density;

    wave_speed = std::sqrt(E/rho);
    dt_max = 0.1 * dx / wave_speed;

    particle_volume = 1.0 / Np; // INITIAL particle volume V^0
    particle_mass = rho * particle_volume;

    particles = Particles(Np);

}




void Simulation::simulate(){

    if (mkdir("dumps", 0777) == -1)
        std::cerr << "Error :  " << strerror(errno) << std::endl;
    else
        std::cout << "Directory created" << std::endl;

    if (mkdir(("dumps/" + sim_name).c_str(), 0777) == -1)
        std::cerr << "Error :  " << strerror(errno) << std::endl;
    else
        std::cout << "Directory created" << std::endl;

    std::ofstream infoFile("dumps/" + sim_name + "/info.txt");
    infoFile << end_frame      << "\n"  // 0
             << frame_dt       << "\n"  // 1
             << dx             << "\n"  // 2
             << mu             << "\n"  // 3
             << lambda         << "\n"  // 4
             << elastic_model  << "\n"  // 5
             << plastic_model  << "\n"  // 6
             << yield_stress   << "\n"  // 7
             << reg_length     << "\n"  // 8
             << reg_const      << "\n";  // 9
    infoFile.close();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Lagrangian coordinates
    particles.x0 = particles.x;
    particles.y0 = particles.y;

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
    // calculateMassConservation();
    explicitEulerUpdate();
    // addExternalParticleGravity();
    G2P();
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

    T max_speed = std::sqrt((particles.vx.array().square() + particles.vy.array().square()).maxCoeff());
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
    T min_x = particles.x.minCoeff();
    T min_y = particles.y.minCoeff();
    T max_x = particles.x.maxCoeff();
    T max_y = particles.y.maxCoeff();

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;

    // ACTUAL midpoints of particles
    T mid_x = min_x + Lx / 2.0;
    T mid_y = min_y + Ly / 2.0;

    // NEW number of grid-dx's per side length
    unsigned int safety_factor = 3;
    Nx = std::ceil(Lx / dx) + safety_factor;
    Ny = std::ceil(Ly / dx) + safety_factor;

    T low_x  = mid_x - Nx*dx/2.0;
    T low_y  = mid_y - Ny*dx/2.0;
    T high_x = mid_x + Nx*dx/2.0;
    T high_y = mid_y + Ny*dx/2.0;

    // NEW grid dimensions N = (L / dx + 1)
    Nx++;
    Ny++;

    debug("               grid   = (", Nx, ", ", Ny, ")");

    // Linspace x and y
    // LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = TVX::LinSpaced(Nx, low_x, high_x);
    grid.y = TVX::LinSpaced(Ny, low_y, high_y);

    grid.vx   = TMX::Zero(Nx, Ny);
    grid.vy   = TMX::Zero(Nx, Ny);
    grid.mass = TMX::Zero(Nx, Ny);

}




void Simulation::moveObjects(T delta_t){
    for (InfinitePlate &obj : objects) {
        obj.move(delta_t);
    }
}

void Simulation::boundaryCollision(T xi, T yi, T& vxi, T& vyi){

    for (InfinitePlate &obj : objects) {
        bool colliding = obj.inside(xi, yi);
        if (colliding) {
            T vx_rel = vxi - obj.vx_object;
            T vy_rel = vyi - obj.vy_object;
            if (boundary_condition == STICKY) { // STICKY
                vx_rel = 0;
                vy_rel = 0;
            } // end STICKY
            vxi = vx_rel + obj.vx_object;
            vyi = vy_rel + obj.vy_object;
        } // end if colliding

    } // end iterator over objects

} // end boundaryCollision


// Currently not working!!
void Simulation::boundaryCorrection(T xi, T yi, T& vxi, T& vyi){

    // trial step
    T x_next = xi + vxi * dt;
    T y_next = yi + vyi * dt;
    moveObjects(dt);

    for (InfinitePlate &obj : objects) {
        bool colliding = obj.inside(x_next, y_next);
        if (colliding) {
            if (obj.plate_type == upper ){
                debug(obj.name, " dist = ", obj.distance(x_next, y_next));
                vyi += obj.distance(x_next, y_next) / dt; // distance is negative since grid point is inside object
            }
            else if (obj.plate_type == lower){
                debug(obj.name, " dist = ", obj.distance(x_next, y_next));
                vyi -= obj.distance(x_next, y_next) / dt;
            }
        } // end if colliding

    } // end iterator over objects

    // Correct back
    moveObjects(-dt);

} // end boundaryCorrection





void Simulation::positionUpdate(){
    for(int p=0; p<Np; p++){
        particles.x(p) = particles.x(p) + dt * particles.vx(p);
        particles.y(p) = particles.y(p) + dt * particles.vy(p);
    } // end loop over particles
}





void Simulation::saveParticleData(std::string extra){
    std::ofstream outFile("dumps/" + sim_name + "/out_part_frame_" + extra + std::to_string(frame) + ".csv");

    outFile << "x"            << ","   // 0
            << "y"            << ","   // 1
            << "z"            << ","   // 2
            << "vx"           << ","   // 3
            << "vy"           << ","   // 4
            << "vz"           << ","   // 5
            << "eps_pl_dev"   << ","   // 6
            << "reg"          << ","   // 7
            << "pressure"     << ","   // 8
            << "devstress"    << "\n"; // 9

    for(int p = 0; p < Np; p++){

        TM2 tau = particles.tau[p];
        T pressure  = -tau.sum() / 2.0;
        TM2 tau_dev = tau + pressure * TM2::Identity();
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        T reg = particles.regularization(p);

        outFile << particles.x(p)              << ","   // 0
                << particles.y(p)              << ","   // 1
                << 0                           << ","   // 2
                << particles.vx(p)             << ","   // 3
                << particles.vy(p)             << ","   // 4
                << 0                           << ","   // 5
                << particles.eps_pl_dev(p)     << ","   // 6
                << reg                         << ","   // 7
                << pressure                    << ","   // 8
                << devstress                   << "\n"; // 9
    }
}

void Simulation::saveGridData(std::string extra){
    std::ofstream outFile("dumps/" + sim_name + "/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    outFile         << "x"       << "," << "y"       << "," << "z" << "," << "vx"         << "," << "vy"         << "," << "vz" << "," << "mass"         << "," << "reg"                    << "\n";
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            outFile << grid.x(i) << "," << grid.y(j) << "," << 0   << "," << grid.vx(i,j) << "," << grid.vy(i,j) << "," << 0    << "," << grid.mass(i,j) << "," << grid.regularization(i,j) << "\n";
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA FUNCTIONS //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function can only be used with fixed mesh
std::pair<TMX, TMX> Simulation::createExternalGridGravity(){
    TMX grid_X0 = TMX::Zero(Nx,Ny);
    TMX grid_Y0 = TMX::Zero(Nx,Ny);
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if (grid.mass(i,j) < 1e-25){ // if no mass at current grid point, no point adding a external force to it.
                grid_X0(i,j) = 0.0;
                grid_Y0(i,j) = 0.0;
            } else {
                grid_X0(i,j) = grid.x(i);
                grid_Y0(i,j) = grid.y(j);
            } // end if grid mass nonzero
        } // end for j
    } // end for i
    return std::make_pair(-2.0*amplitude*grid_X0, -2.0*amplitude*grid_Y0);
} // end createExternalGridGravity


// Much more computationally expensive than createExternalGridGravity,
// but can be used with remeshing
void Simulation::addExternalParticleGravity(){
    // 1. Transfer grid velocity to particles
    G2P();
    // 2. Apply gravity on particle velocity
    for(int p=0; p<Np; p++){
        particles.vx(p) += dt * (-2.0*amplitude*particles.x0(p));
        particles.vy(p) += dt * (-2.0*amplitude*particles.y0(p));
    }
    // 3. Transfer particle velocity back to grid
    P2G();
} // end addExternalParticleGravity


// These functions are only for validating conservation laws
void Simulation::calculateMomentumOnParticles(){
    T momentum_x = 0;
    T momentum_y = 0;
    for(int p=0; p<Np; p++){
        momentum_x += particle_mass * particles.vx(p);
        momentum_y += particle_mass * particles.vy(p);
    }
    debug("               total part momentum x-comp = ", momentum_x);
    debug("               total part momentum y-comp = ", momentum_y);
}
void Simulation::calculateMomentumOnGrid(){
    T momentum_x = 0;
    T momentum_y = 0;
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            momentum_x += grid.mass(i,j) * grid.vx(i,j);
            momentum_y += grid.mass(i,j) * grid.vy(i,j);
        }
    }
    debug("               total grid momentum x-comp = ", momentum_x);
    debug("               total grid momentum y-comp = ", momentum_y);
}

void Simulation::calculateMassConservation(){
    T grid_mass_total     = grid.mass.sum();
    T particle_mass_total = particle_mass*Np;
    debug("               total grid mass = ", grid_mass_total    );
    debug("               total part mass = ", particle_mass_total);
    if ( std::abs(grid_mass_total-particle_mass_total) > 1e-5 * particle_mass_total ){
        debug("MASS NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}
