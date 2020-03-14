#include "simulation.hpp"
#include <omp.h>
int N_THREADS = 4;

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

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Lagrangian coordinates
    particles.x0 = particles.x;
    particles.y0 = particles.y;

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;
    saveSim();
    saveGridVelocities();
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
            saveSim();
            saveGridVelocities();
        }
        if (std::abs(final_time-time) < 1e-15 || final_time < time){
            std::cout << "The simulation ended successfully at time = " << time << std::endl;
            break;
        }
    }
    saveSim();
    saveGridVelocities();

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
    moveObjects();
    P2G();
    // calculateMassConservation();
    explicitEulerUpdate();
    // addExternalParticleGravity();
    G2P();
    deformationUpdate();
    positionUpdate();
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
    explicitEulerUpdate_Optimized();
    t_euler.stop(); runtime_euler += t_euler.get_timing();
}

void Simulation::G2P(){
    timer t_g2p; t_g2p.start();
    // G2P_Baseline();
    // G2P_Optimized();
    G2P_Optimized_Parallel();
    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();
}


void Simulation::P2G_Baseline(){
    // This (nested) loop over i (and j) can easily be paralellized
    // Need to create thread-local grid.mass, grid.vx, grid.vy
    for(int i=0; i<Nx; i++){
        T xi = grid.x(i);
        for(int j=0; j<Ny; j++){
            T yi = grid.y(j);
            T vxi = 0;
            T vyi = 0;
            T mass = 0;
            for(int p=0; p<Np; p++){
                T xp = particles.x(p);
                T yp = particles.y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, dx);
                    mass += weight;
                    vxi  += particles.vx(p) * weight;
                    vyi  += particles.vy(p) * weight;
                }
            } // end for particles
            grid.mass(i,j) = mass * particle_mass;
            if (mass < 1e-25){
                grid.vx(i,j)   = 0.0;
                grid.vy(i,j)   = 0.0;
            } else {
                grid.vx(i,j)   = vxi / mass;
                grid.vy(i,j)   = vyi / mass;
            }
        } // end for j
    } // end for i
} // end P2G_Baseline


void Simulation::P2G_Optimized(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    for(int p = 0; p < Np; p++){
        T xp = particles.x(p);
        T yp = particles.y(p);
        unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);
                T weight = wip(xp, yp, xi, yi, dx);

                if (weight > 1e-25){
                    grid.mass(i,j) += weight;
                    grid.vx(i,j)   += particles.vx(p) * weight;
                    grid.vy(i,j)   += particles.vy(p) * weight;
                }

            } // end for j
        } // end for i
    } // end for p

    ///////////////////////////////////////////////////////////
    // At this point in time grid.mass is equal to m_i / m_p //
    ///////////////////////////////////////////////////////////
    grid.vx = (grid.mass.array() > 0).select( (grid.vx.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.vy = (grid.mass.array() > 0).select( (grid.vy.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.mass *= particle_mass;

} // end P2G_Optimized


void Simulation::P2G_Optimized_Parallel(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    grid.mass.setZero(Nx, Ny);
    grid.vx.setZero(Nx, Ny);
    grid.vy.setZero(Nx, Ny);

    #pragma omp parallel num_threads(N_THREADS)
    {
        TMX grid_mass_local = TMX::Zero(Nx, Ny);
        TMX grid_vx_local   = TMX::Zero(Nx, Ny);
        TMX grid_vy_local   = TMX::Zero(Nx, Ny);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            T xp = particles.x(p);
            T yp = particles.y(p);
            unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x(i);
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y(j);
                    T weight = wip(xp, yp, xi, yi, dx);

                    if (weight > 1e-25){
                        grid_mass_local(i,j) += weight;
                        grid_vx_local(i,j)   += particles.vx(p) * weight;
                        grid_vy_local(i,j)   += particles.vy(p) * weight;
                    }

                } // end for j
            } // end for i
        } // end for p

        #pragma omp critical
        {
            for(int i = 0; i < Nx; i++){
                for(int j = 0; j < Ny; j++){
                    grid.mass(i,j) += grid_mass_local(i,j);
                    grid.vx(i,j)   += grid_vx_local(i,j);
                    grid.vy(i,j)   += grid_vy_local(i,j);
                } // end for j
            } // end for i
        } // end omp critical

    } // end omp parallel

    ///////////////////////////////////////////////////////////
    // At this point in time grid.mass is equal to m_i / m_p //
    ///////////////////////////////////////////////////////////
    grid.vx = (grid.mass.array() > 0).select( (grid.vx.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.vy = (grid.mass.array() > 0).select( (grid.vy.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.mass *= particle_mass;

} // end P2G_Optimized_Parallel






void Simulation::explicitEulerUpdate_Baseline(){
    TV2 grad_wip;
    TM2 Fe, dPsidF;

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    // TV2 external_gravity;
    ////////////////////////////////////////////////////////

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if (grid.mass(i,j) > 1e-25){
                T xi = grid.x(i);
                T yi = grid.y(j);
                TV2 grid_force = TV2::Zero();
                for(int p=0; p<Np; p++){
                    T xp = particles.x(p);
                    T yp = particles.y(p);
                    if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                        // Fe = particles[p].F;
                        Fe = particles.F[p];

                        // Remember:
                        // P = dPsidF              (first Piola-Kirchhoff stress tensor)
                        // tau = P * F.transpose() (Kirchhoff stress tensor)

                        if (elastic_model == NeoHookean){
                            dPsidF = mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
                        }
                        else if (elastic_model == StvkWithHencky){
                            Eigen::JacobiSVD<TM2> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
                            TA2 sigma = svd.singularValues().array(); // abs() for inverse also??
                            TM2 logSigma = sigma.abs().log().matrix().asDiagonal();
                            TM2 invSigma = sigma.inverse().matrix().asDiagonal();
                            dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
                        }
                        else{
                            debug("You specified an unvalid ELASTIC model!");
                        }

                        grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                        grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                        grid_force += dPsidF * Fe.transpose() * grad_wip;

                    }
                } // end for particles

                TV2 velocity_increment = -dt * particle_volume * grid_force / grid.mass(i,j) + dt * gravity;

                //////////// if external grid gravity: //////////////////
                // external_gravity(0) = external_gravity_pair.first(i,j);
                // external_gravity(1) = external_gravity_pair.second(i,j);
                // velocity_increment += dt * external_gravity;
                ////////////////////////////////////////////////////////

                T new_vxi = grid.vx(i,j) + velocity_increment(0);
                T new_vyi = grid.vy(i,j) + velocity_increment(1);
                T new_xi = grid.x(i) + dt * new_vxi;
                T new_yi = grid.y(j) + dt * new_vyi;
                boundaryCollision(new_xi, new_yi, new_vxi, new_vyi);

                grid.vx(i,j) = new_vxi;
                grid.vy(i,j) = new_vyi;
            } // end if positive mass
        } // end for j
    } // end for i
} // end explicitEulerUpdate_Baseline

void Simulation::explicitEulerUpdate_Optimized(){
    TV2 grad_wip;
    TM2 Fe, dPsidF, tau;

    // Remember:
    // P = dPsidF              (first Piola-Kirchhoff stress tensor)
    // tau = P * F.transpose() (Kirchhoff stress tensor)

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    TMX grid_force_x = TMX::Zero(Nx, Ny);
    TMX grid_force_y = TMX::Zero(Nx, Ny);

    for(int p = 0; p < Np; p++){

        Fe = particles.F[p];

        if (elastic_model == NeoHookean){
            dPsidF = mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
        }
        else if (elastic_model == StvkWithHencky){ // St Venant Kirchhoff with Hencky strain
            Eigen::JacobiSVD<TM2> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
            TA2 sigma = svd.singularValues().array(); // abs() for inverse also??
            TM2 logSigma = sigma.abs().log().matrix().asDiagonal();
            TM2 invSigma = sigma.inverse().matrix().asDiagonal();
            dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
        }
        else{
            debug("You specified an unvalid ELASTIC model!");
        }

        tau = dPsidF * Fe.transpose();
        particles.tau[p] = tau;

        T xp = particles.x(p);
        T yp = particles.y(p);
        unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);

                if (grid.mass(i,j) > 1e-25){
                    grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                    grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);
                    TV2 grid_force_increment = tau * grad_wip;
                    grid_force_x(i,j) += grid_force_increment(0);
                    grid_force_y(i,j) += grid_force_increment(1);
                } // end if non-zero grid mass

             } // end for j
         } // end for i

    } // end for particles

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    ////////////////////////////////////////////////////////

    T dt_particle_volume = dt * particle_volume;
    T dt_gravity_x = dt * gravity(0);
    T dt_gravity_y = dt * gravity(1);

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            T mi = grid.mass(i,j);
            if (mi > 1e-25){

                T velocity_increment_x = -dt_particle_volume * grid_force_x(i,j) / mi + dt_gravity_x;
                T velocity_increment_y = -dt_particle_volume * grid_force_y(i,j) / mi + dt_gravity_y;

                //////////// if external grid gravity: //////////////////
                // T external_gravity = external_gravity_pair.first(i,j);
                // T external_gravity = external_gravity_pair.second(i,j);
                // velocity_increment_x += dt * external_gravity(0);
                // velocity_increment_y += dt * external_gravity(1);
                ////////////////////////////////////////////////////////

                T new_vxi = grid.vx(i,j) + velocity_increment_x;
                T new_vyi = grid.vy(i,j) + velocity_increment_y;
                T new_xi = grid.x(i) + dt * new_vxi;
                T new_yi = grid.y(j) + dt * new_vyi;
                boundaryCollision(new_xi, new_yi, new_vxi, new_vyi);
                grid.vx(i,j) = new_vxi;
                grid.vy(i,j) = new_vyi;

            } // end if non-zero grid mass

        } // end for j
    } // end for i

} // end explicitEulerUpdate_Optimized


void Simulation::moveObjects(){
    for (InfinitePlate &obj : objects) {
        obj.move(dt);
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





void Simulation::G2P_Baseline(){
    // This loop over p can be easily paralellized
    // Need to create thread-local versions of particles.vx and vy
    for(int p=0; p<Np; p++){
        T xp = particles.x(p);
        T yp = particles.y(p);
        T vxp = 0;
        T vyp = 0;
        for(int i=0; i<Nx; i++){
            T xi = grid.x(i);
            for(int j=0; j<Ny; j++){
                T yi = grid.y(j);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid.vx(i,j) * weight;
                    vyp += grid.vy(i,j) * weight;
                }
            }
        }
        particles.vx(p) = vxp;
        particles.vy(p) = vyp;
    }
} // end G2P_Baseline

void Simulation::G2P_Optimized(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    // This loop over p can be easily paralellized
    // Need to create thread-local versions of particles.vx and vy
    for(int p = 0; p < Np; p++){
        T xp = particles.x(p);
        T yp = particles.y(p);
        T vxp = 0;
        T vyp = 0;
        unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);
                T weight = wip(xp, yp, xi, yi, dx);
                vxp += grid.vx(i,j) * weight;
                vyp += grid.vy(i,j) * weight;
            } // end loop j
        } // end loop i
        particles.vx(p) = vxp;
        particles.vy(p) = vyp;
    } // end loop p
} // end G2P_Optimized

void Simulation::G2P_Optimized_Parallel(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    particles.vx.setZero(Np);
    particles.vy.setZero(Np);

    #pragma omp parallel num_threads(N_THREADS)
    {
        TVX particles_vx_local = TVX::Zero(Np);
        TVX particles_vy_local = TVX::Zero(Np);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            T xp = particles.x(p);
            T yp = particles.y(p);
            T vxp = 0;
            T vyp = 0;
            unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x(i);
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y(j);
                    T weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid.vx(i,j) * weight;
                    vyp += grid.vy(i,j) * weight;
                } // end loop j
            } // end loop i
            particles_vx_local(p) = vxp;
            particles_vy_local(p) = vyp;
        } // end loop p

        #pragma omp critical
        {
            for(int p = 0; p < Np; p++){
                particles.vx(p) += particles_vx_local(p);
                particles.vy(p) += particles_vy_local(p);
            }
        } // end omp critical

    } // end omp paralell

} // end G2P_Optimized_Parallel






// Deformation gradient is updated based on the NEW GRID VELOCITIES and the OLD PARTICLE POSITIONS
void Simulation::deformationUpdate(){
    timer t_defgrad; t_defgrad.start();

    unsigned int plastic_count = 0;

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    for(int p=0; p<Np; p++){

        TM2 sum = TM2::Zero();

        T xp = particles.x(p);
        T yp = particles.y(p);

        unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);

                TV2 vi;
                vi(0) = grid.vx(i,j);
                vi(1) = grid.vy(i,j);

                TV2 grad_wip;
                grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                sum += vi * grad_wip.transpose();
            } // end loop i
        } // end loop j

        TM2 Fe_trial = particles.F[p];
        Fe_trial = Fe_trial + dt * sum * Fe_trial;
        particles.F[p] = Fe_trial;

        if (plastic_model == VonMises){
            Eigen::JacobiSVD<TM2> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
            TV2 hencky = svd.singularValues().array().log(); // Jixie does not use abs value, however Pradhana-thesis does.
            T   hencky_trace = hencky.sum();
            TV2 hencky_deviatoric = hencky - (hencky_trace / 2.0) * TV2::Ones();
            T   hencky_deviatoric_norm = hencky_deviatoric.norm();

            T delta_gamma = hencky_deviatoric_norm - yield_stress / (2 * mu);
            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.eps_pl_dev[p] += delta_gamma;
                hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        } // end VonMises Plasticity
        else if (plastic_model == NoPlasticity){
            // Do nothing
        }
        else{
            debug("You specified an unvalid PLASTIC model!");
        }
    } // end loop over particles

    debug("               projected particles = ", plastic_count, " / ", Np);
    t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();

} // end deformationUpdate





void Simulation::positionUpdate(){
    for(int p=0; p<Np; p++){
        particles.x(p) = particles.x(p) + dt * particles.vx(p);
        particles.y(p) = particles.y(p) + dt * particles.vy(p);
    } // end loop over particles
}





void Simulation::saveSim(std::string extra){
    std::ofstream outFile("dumps/out_part_frame_" + extra + std::to_string(frame) + ".csv");
    for(int p = 0; p < Np; p++){

        TM2 tau = particles.tau[p];
        T pressure  = -tau.sum() / 2.0;
        TM2 tau_dev = tau + pressure * TM2::Identity();
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        outFile << particles.x[p]          << ","
                << particles.y[p]          << ","
                << 0                       << ","
                << particles.vx[p]         << ","
                << particles.vy[p]         << ","
                << 0                       << ","
                << particles.eps_pl_dev[p] << ","
                << pressure                << ","
                << devstress               << "\n";
    }
}

void Simulation::saveGridVelocities(std::string extra){
    std::ofstream outFile("dumps/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            outFile << grid.x(i) << "," << grid.y(j) << "," << 0 << "," << grid.vx(i,j) << "," << grid.vy(i,j) << "," << 0 << "\n";
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////??????/ EXTRA FUNCTIONS //////////////////////////////////////////////////
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
