#include "simulation.hpp"
#include "tools.hpp"

void Simulation::setElasticParams(double E, double nu, double density){
    mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
    lambda = E / (2.0*(1.0+nu));
    rho = density;
    dt_max = 0.5 * dx * std::sqrt(rho/E);
}

void Simulation::simulate(){
    double t = 0;
    for (int i = 0; i < Nt; i++){
        advanceStep();
        t += dt;
        current_step++;
        std::cout << "Step: " << current_step << "\t Time: " << t << std::endl;
        saveSim();
    }
}

void Simulation::advanceStep(){
    remesh();
    P2G();
    explicitEulerUpdate();
    G2P();
    deformationUpdate();
    positionUpdate();
}




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

    //std::cout << grid_X.rows() << std::endl;
    //std::cout << grid_X.cols() << std::endl;
    grid_VX   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_VY   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_mass = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());

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
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
                    double weight = wip(xp, yp, xi, yi, dx);
                    mass += weight;
                    vxi  += particles_vx(p) * weight;
                    vyi  += particles_vy(p) * weight;
                }
            }
            //debug(vxi);
            grid_mass(i,j) = mass * particle_mass;
            grid_VX(i,j)   = vxi  * particle_mass / mass;
            grid_VY(i,j)   = vyi  * particle_mass / mass;
        }
    }
}





void Simulation::explicitEulerUpdate(){
    TV2 grad_wip;
    Eigen::Array2d sigma;
    TM2 Fe, logSigma, invSigma, dPsidF;
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            double xi = grid_X(i,j);
            double yi = grid_Y(i,j);
            TV2 force = TV2::Zero();
            for(int p=0; p<Np; p++){
                double xp = particles_x(p);
                double yp = particles_y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
                    Fe = particles[p].F; // must update F!!

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
            }
            TV2 velocity_increment = dt * force / grid_mass(i,j);
            grid_VX(i,j) += velocity_increment(0);
            grid_VY(i,j) += velocity_increment(1);
        }
    }
}





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
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
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
    } // end loop over particles


}






void Simulation::saveSim(){
  std::ofstream outFile("out_" + std::to_string(current_step) + ".csv");
  for(int p = 0; p < Np; p++){
      outFile << p << "," << particles_x[p] << "\t" << particles_y[p] << "\t" << particles_vx[p] << "\t" << particles_vy[p] << "\n";
  }
}





Simulation::Simulation() : current_step(0) {
  // default case: unit box with 10 times dx = 0.1
    Nt = 2;
    dx = 0.1;
    rho = 700;
    setElasticParams(1e5, 0.3, rho);
    debug("dt_max = ", dt_max);
    dt = 0.002;

    Nx = 10;
    Ny = 10;
    Np = Nx*Ny*4;
    debug("Np = ", Np);

    particle_volume = 1.0 / Np;
    particle_mass = rho * 1.0 / Np;

    grid_X = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_Y = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_VX = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_VY = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_mass = Eigen::MatrixXd::Zero(Nx, Ny);

    Particle part;
    particles.resize(Np);
    std::fill(particles.begin(), particles.end(), part);

    particles_x  = Eigen::VectorXd::Zero(Np);
    particles_y  = Eigen::VectorXd::Zero(Np);
    particles_vx = Eigen::VectorXd::Zero(Np);
    particles_vy = Eigen::VectorXd::Zero(Np);

    double amplitude = 0.01;

    int p = -1;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){

            p++;
            double px = (i+0.25)*dx;
            double py = (j+0.25)*dx;
            double pvx = amplitude*std::sin(M_PI*px);
            double pvy = amplitude*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.75)*dx;
            py = (j+0.75)*dx;
            pvx = amplitude*std::sin(M_PI*px);
            pvy = amplitude*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.25)*dx;
            py = (j+0.75)*dx;
            pvx = amplitude*std::sin(M_PI*px);
            pvy = amplitude*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.75)*dx;
            py = (j+0.25)*dx;
            pvx = amplitude*std::sin(M_PI*px);
            pvy = amplitude*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;
        } // end for i
    } // end for j
}
