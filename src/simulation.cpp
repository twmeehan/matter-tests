#include "simulation.hpp"

Simulation::Simulation(){
    current_time_step = 0;
    time = 0;
    frame = 0;
    dt = 1e-3; // just initial guess
    exit = 0;

    // N = (L / dx + 1)
    Nx = 21;
    Ny = 21;

    ////////////// If remesh every time step /////////////////
    grid_X    = TMX::Zero(Nx, Ny);
    grid_Y    = TMX::Zero(Nx, Ny);
    grid_VX   = TMX::Zero(Nx, Ny);
    grid_VY   = TMX::Zero(Nx, Ny);
    grid_mass = TMX::Zero(Nx, Ny);
    ///////////////// If constant mesh ////////////////////////
    /*
    TVX lin_x = TVX::LinSpaced(sim.Nx, -0.5, 1.5);
    TVX lin_y = TVX::LinSpaced(sim.Ny, -0.5, 1.5);

    sim.grid_X.resize(sim.Ny, sim.Nx);
    sim.grid_Y.resize(sim.Ny, sim.Nx);
    for (int i = 0; i < sim.Ny; ++i) {
      sim.grid_X.row(i) = lin_x.transpose();
    }
    for (int j = 0; j < sim.Nx; ++j) {
      sim.grid_Y.col(j) = lin_y;
    }
    sim.grid_X.transposeInPlace();
    sim.grid_Y.transposeInPlace();

    sim.grid_VX   = TMX::Zero(sim.Nx, sim.Ny);
    sim.grid_VY   = TMX::Zero(sim.Nx, sim.Ny);
    sim.grid_mass = TMX::Zero(sim.Nx, sim.Ny);
    */
    ////////////////////////////////////////////////////////////////
}

void Simulation::initialize(T E, T nu, T density){
    mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    lambda = E / (2.0*(1.0+nu));
    rho = density;

    wave_speed = std::sqrt(E/rho);
    dt_max = 0.1 * dx / wave_speed ;

    particle_volume = 1.0 / Np; // INITIAL particle volume V^0
    particle_mass = rho * particle_volume;

    // Particle part;
    // sim.particles.resize(sim.Np);
    // std::fill(sim.particles.begin(), sim.particles.end(), part);
    particles_F.resize(Np);
    std::fill(particles_F.begin(), particles_F.end(), TM2::Identity());

    particles_x  = TVX::Zero(Np);
    particles_y  = TVX::Zero(Np);
    particles_vx = TVX::Zero(Np);
    particles_vy = TVX::Zero(Np);
}

void Simulation::simulate(){
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
}

void Simulation::advanceStep(){
    updateDt();
    remesh();
    P2G();
    explicitEulerUpdate();
    G2P();
    deformationUpdate();
    positionUpdate();
}




void Simulation::updateDt(){

    T max_speed = std::sqrt((particles_vx.array().square() + particles_vy.array().square()).maxCoeff());
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
    T min_x = particles_x.minCoeff();
    T min_y = particles_y.minCoeff();
    T max_x = particles_x.maxCoeff();
    T max_y = particles_y.maxCoeff();

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;

    // ACTUAL midpoints of particles
    T mid_x = min_x + Lx / 2.0;
    T mid_y = min_y + Ly / 2.0;

    // NEW number of grid-dx's per side length
    unsigned int safety_factor = 1;
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
    TVX lin_x = TVX::LinSpaced(Nx, low_x, high_x);
    TVX lin_y = TVX::LinSpaced(Ny, low_y, high_y);

    // Meshgrid
    grid_X.resize(Ny, Nx);
    grid_Y.resize(Ny, Nx);
    for (int i = 0; i < Ny; ++i) {
      grid_X.row(i) = lin_x.transpose();
    }
    for (int j = 0; j < Nx; ++j) {
      grid_Y.col(j) = lin_y;
    }
    grid_X.transposeInPlace();
    grid_Y.transposeInPlace();

    grid_VX   = TMX::Zero(grid_X.rows(), grid_X.cols());
    grid_VY   = TMX::Zero(grid_X.rows(), grid_X.cols());
    grid_mass = TMX::Zero(grid_X.rows(), grid_X.cols());

}




void Simulation::P2G(){
    for(int i=0; i<Nx; i++){
        //debug("i = ", i);
        for(int j=0; j<Ny; j++){
            //debug("j = ", j);
            T xi = grid_X(i,j);
            T yi = grid_Y(i,j);
            T vxi = 0;
            T vyi = 0;
            T mass = 0;
            for(int p=0; p<Np; p++){
                //debug(p);
                T xp = particles_x(p);
                T yp = particles_y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, dx);
                    mass += weight;
                    vxi  += particles_vx(p) * weight;
                    vyi  += particles_vy(p) * weight;
                }
            } // end for particles
            grid_mass(i,j) = mass * particle_mass;
            if (mass < 1e-15){
                grid_VX(i,j)   = 0.0;
                grid_VY(i,j)   = 0.0;
            } else {
                grid_VX(i,j)   = vxi / mass;
                grid_VY(i,j)   = vyi / mass;
            }
        } // end for j
    } // end for i

    // if (current_time_step == 0)
    //     debug(grid_VX);

    debug("               total grid mass = ", grid_mass.sum());
    debug("               total part mass = ", particle_mass*Np);


} // end P2G





void Simulation::explicitEulerUpdate(){
    TV2 grad_wip;
    TM2 Fe, dPsidF;

    // Only for St Venant Kirchhoff
    TM2 logSigma, invSigma; TA2 sigma;

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if (grid_mass(i,j) > 1e-15){
                T xi = grid_X(i,j);
                T yi = grid_Y(i,j);
                TV2 force = TV2::Zero();
                for(int p=0; p<Np; p++){
                    T xp = particles_x(p);
                    T yp = particles_y(p);
                    if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                        // Fe = particles[p].F;
                        Fe = particles_F[p];

                        // Remember:
                        // P = dPsidF              (first Piola-Kirchhoff stress tensor)
                        // tau = P * F.transpose() (Kirchhoff stress tensor)

                        /////// St Venant Kirchhoff with Hencky strain ///////////////
                        Eigen::JacobiSVD<TM2> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
                        sigma = svd.singularValues().array(); // abs() for inverse also??
                        logSigma = sigma.abs().log().matrix().asDiagonal();
                        invSigma = sigma.inverse().matrix().asDiagonal();
                        dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
                        ///////////////////////////////////////////////////////////////

                        ////////////// Neo-Hookean ////////////////////////////////////
                        // dPsidF = mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
                        ///////////////////////////////////////////////////////////////

                        grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                        grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                        force += dPsidF * Fe.transpose() * grad_wip;
                    }
                } // end for particles
                TV2 velocity_increment = -dt * particle_volume * force / grid_mass(i,j) + dt * gravity;
                grid_VX(i,j) += velocity_increment(0);
                grid_VY(i,j) += velocity_increment(1);
            } // end if positive mass
        } // end for j
    } // end for i
} // end explicitEulerUpdate





void Simulation::G2P(){
    for(int p=0; p<Np; p++){
        T xp = particles_x(p);
        T yp = particles_y(p);
        T vxp = 0;
        T vyp = 0;
        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){
                T xi = grid_X(i,j);
                T yi = grid_Y(i,j);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid_VX(i,j) * weight;
                    vyp += grid_VY(i,j) * weight;
                }
            }
        }
        particles_vx(p) = vxp;
        particles_vy(p) = vyp;
    }
}






void Simulation::deformationUpdate(){
    unsigned int plastic_count = 0;
    for(int p=0; p<Np; p++){

        TM2 sum = TM2::Zero();

        T xp = particles_x(p);
        T yp = particles_y(p);

        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){

                T xi = grid_X(i,j);
                T yi = grid_Y(i,j);

                TV2 vi;
                vi(0) = grid_VX(i,j);
                vi(1) = grid_VY(i,j);

                TV2 grad_wip;
                grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                sum += vi * grad_wip.transpose();
            } // end loop i
        } // end loop j

        // TM2 Fe = particles[p].F;
        // particles[p].F = Fe + dt * sum * Fe;
        TM2 Fe_trial = particles_F[p];
        Fe_trial = Fe_trial + dt * sum * Fe_trial;
        particles_F[p] = Fe_trial;

        if (plasticity){
            Eigen::JacobiSVD<TM2> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
            TV2 hencky_trial = svd.singularValues().array().log(); // Jixie does not use abs value, however Pradhana-thesis does.
            T   hencky_trial_trace = hencky_trial.sum();
            TV2 hencky_trial_deviatoric = hencky_trial - (hencky_trial_trace / 2.0) * TV2::Ones();
            T   hencky_trial_deviatoric_norm = hencky_trial_deviatoric.norm();

            // Von Mises:
            T delta_gamma = hencky_trial_deviatoric_norm - yield_stress / (2 * mu);
            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                TV2 hencky_new = hencky_trial - delta_gamma * (hencky_trial_deviatoric / hencky_trial_deviatoric_norm);
                particles_F[p] = svd.matrixU() * hencky_new.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        }

    } // end loop over particles

    debug("               projected particles = ", plastic_count, " / ", Np);

} // end deformationUpdate





void Simulation::positionUpdate(){
    for(int p=0; p<Np; p++){
        particles_x(p) = particles_x(p) + dt * particles_vx(p);
        particles_y(p) = particles_y(p) + dt * particles_vy(p);

        if (particles_x(p) > 1.5 || particles_x(p) < -0.5 || particles_y(p) > 1.5 || particles_y(p) < -0.5){
            debug("Particle goes outside bounding box!!!");
        }
    } // end loop over particles


}






void Simulation::saveSim(std::string extra){
    std::ofstream outFile("dumps/out_part_frame_" + extra + std::to_string(frame) + ".csv");
    for(int p = 0; p < Np; p++){
        outFile << particles_x[p] << "," << particles_y[p] << "," << 0 << "," << particles_vx[p] << "," << particles_vy[p] << "," << 0 << "\n";
    }
}

void Simulation::saveGridVelocities(std::string extra){
    std::ofstream outFile("dumps/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            int k = i*Ny+j;
            outFile << grid_X(i,j) << "," << grid_Y(i,j) << "," << 0 << "," << grid_VX(i,j) << "," << grid_VY(i,j) << "," << 0 << "\n";
        }
    }
}
