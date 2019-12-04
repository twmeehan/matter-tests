#include "simulation.hpp"
#include "tools.hpp"

void Simulation::setElasticParams(double E, double nu, double density){
    mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    lambda = E / (2.0*(1.0+nu));
    rho = density;
    dt_max = 0.5 * dx * std::sqrt(rho/E);
}

void Simulation::simulate(){
    saveSim();
    saveGridVelocities();
    double t = 0;
    for (int i = 0; i < max_time_steps; i++){
        advanceStep();
        if (exit == 1)
            return;
        t += dt;
        current_time_step++;
        std::cout << "Step: " << current_time_step << "    Time: " << t << std::endl;
        saveSim();
        saveGridVelocities();
        if (t >= T){
            std::cout << "The simulation ended successfully at time t = " << t << std::endl;
            break;
        }
    }
    saveSim();
    saveGridVelocities();
}

void Simulation::advanceStep(){
    updateDt();
    //remesh();
    P2G();
    explicitEulerUpdate();
    G2P();
    deformationUpdate();
    positionUpdate();
}




void Simulation::updateDt(){

    double max_speed = std::sqrt((particles_vx.array().square() + particles_vy.array().square()).maxCoeff());
    double dt_cfl = cfl * dx / max_speed;
    dt = std::min(dt_cfl, dt_max);
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
} // end updateDt




void Simulation::remesh(){

    double min_x = particles_x.minCoeff();
    double min_y = particles_y.minCoeff();
    double max_x = particles_x.maxCoeff();
    double max_y = particles_y.maxCoeff();

    double Lx = max_x - min_x;
    double Ly = max_y - min_y;

    Nx = std::ceil(Lx / dx);
    Ny = std::ceil(Ly / dx);

    double mid_x = min_x + Lx / 2.0;
    double mid_y = min_y + Ly / 2.0;

    double low_x  = mid_x - Nx*dx/2.0;
    double low_y  = mid_y - Ny*dx/2.0;
    double high_x = mid_x + Nx*dx/2.0;
    double high_y = mid_y + Ny*dx/2.0;

    // N = (L / dx + 1)
    Nx++;
    Ny++;

    debug("               grid   = (", Nx, ", ", Ny, ")");

    // Linspace x and y
    Eigen::VectorXd lin_x = Eigen::VectorXd::LinSpaced(Nx, low_x, high_x);
    Eigen::VectorXd lin_y = Eigen::VectorXd::LinSpaced(Ny, low_y, high_y);

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

    // if (current_time_step==0)
    //     debug(grid_Y);

    grid_VX   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_VY   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_mass = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());

    // debug("hello");

}




void Simulation::P2G(){
    for(int i=0; i<Nx; i++){
        //debug("i = ", i);
        for(int j=0; j<Ny; j++){
            //debug("j = ", j);
            double xi = grid_X(i,j);
            double yi = grid_Y(i,j);
            double vxi = 0;
            double vyi = 0;
            double mass = 0;
            for(int p=0; p<Np; p++){
                //debug(p);
                double xp = particles_x(p);
                double yp = particles_y(p);
                if ( std::abs(xp-xi) < 2*dx || std::abs(xp-yi) < 2*dx){
                    double weight = wip(xp, yp, xi, yi, dx);
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
                grid_VX(i,j)   = vxi  * particle_mass / mass;
                grid_VY(i,j)   = vyi  * particle_mass / mass;
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
    Eigen::Array2d sigma;
    TM2 Fe, logSigma, invSigma, dPsidF;
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if (grid_mass(i,j) > 1e-15){
                double xi = grid_X(i,j);
                double yi = grid_Y(i,j);
                TV2 force = TV2::Zero();
                for(int p=0; p<Np; p++){
                    double xp = particles_x(p);
                    double yp = particles_y(p);
                    if ( std::abs(xp-xi) < 2*dx || std::abs(xp-yi) < 2*dx){
                        Fe = particles[p].F;

                        Eigen::JacobiSVD<TM2> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
                        sigma = svd.singularValues().array().abs().max(1e-4);

                        // should optimize code by working in principal space

                        logSigma = sigma.log().matrix().asDiagonal();
                        invSigma = sigma.inverse().matrix().asDiagonal();

                        // debug("sigma = \n", sigma);
                        // debug("logSigma = \n", logSigma);
                        // debug("invSigma = \n", invSigma);

                        dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();

                        grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                        grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                        force -= particle_volume * dPsidF * Fe.transpose() * grad_wip; // pull out particle volume and add gravity
                    }
                } // end for particles
                TV2 velocity_increment = dt * force / grid_mass(i,j);
                grid_VX(i,j) += velocity_increment(0);
                grid_VY(i,j) += velocity_increment(1);
            } // end if positive mass
        } // end for j
    } // end for i
} // end explicitEulerUpdate





void Simulation::G2P(){
    for(int p=0; p<Np; p++){
        double xp = particles_x(p);
        double yp = particles_y(p);
        double vxp = 0;
        double vyp = 0;
        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){
                double xi = grid_X(i,j);
                double yi = grid_Y(i,j);
                if ( std::abs(xp-xi) < 1.5*dx || std::abs(xp-yi) < 1.5*dx){
                    double weight = wip(xp, yp, xi, yi, dx);
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
    for(int p=0; p<Np; p++){

        TM2 Fe = particles[p].F;
        TM2 Fe_new = Fe;

        double xp = particles_x(p);
        double yp = particles_y(p);

        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){

                double xi = grid_X(i,j);
                double yi = grid_Y(i,j);

                TV2 vi;
                vi(0) = grid_VX(i,j);
                vi(1) = grid_VY(i,j);

                TV2 grad_wip;
                grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);

                Fe_new += dt * (vi * grad_wip.transpose()) * Fe;
            } // end loop i
        } // end loop j

        particles[p].F = Fe_new;

    } // end loop over particles

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






void Simulation::saveSim(){
    std::ofstream outFile("out_" + std::to_string(current_time_step) + ".csv");
    for(int p = 0; p < Np; p++){
        outFile << p << "," << particles_x[p] << "," << particles_y[p] << "," << particles_vx[p] << "," << particles_vy[p] << "\n";
    }
}

void Simulation::saveGridVelocities(){
    std::ofstream outFile("out_gridvel_" + std::to_string(current_time_step) + ".csv");
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            int k = i*Ny+j;
            outFile << k << "," << grid_X(i,j) << "," << grid_Y(i,j) << "," << grid_VX(i,j) << "," << grid_VY(i,j) << "\n";
        }
    }
}
