#include "simulation.hpp"
#include <omp.h>

void Simulation::explicitEulerUpdate_Optimized_Parallel(){
    TV2 grad_wip;
    TM2 Fe, dPsidF, tau;

    // Remember:
    // P = dPsidF              (first Piola-Kirchhoff stress tensor)
    // tau = P * F.transpose() (Kirchhoff stress tensor)

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    TMX grid_force_x = TMX::Zero(Nx, Ny);
    TMX grid_force_y = TMX::Zero(Nx, Ny);

    #pragma omp parallel num_threads(n_threads)
    {
        TMX grid_force_x_local = TMX::Zero(Nx, Ny);
        TMX grid_force_y_local = TMX::Zero(Nx, Ny);

        #pragma omp for
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
                        grid_force_x_local(i,j) += grid_force_increment(0);
                        grid_force_y_local(i,j) += grid_force_increment(1);
                    } // end if non-zero grid mass

                 } // end for j
             } // end for i

        } // end for particles

        #pragma omp critical
        {
            for(int i = 0; i < Nx; i++){
                for(int j = 0; j < Ny; j++){
                    grid_force_x(i,j) += grid_force_x_local(i,j);
                    grid_force_y(i,j) += grid_force_y_local(i,j);
                } // end for j
            } // end for i
        } // end omp critical

    } // end omp parallel

    T dt_particle_volume = dt * particle_volume;
    T dt_gravity_x = dt * gravity(0);
    T dt_gravity_y = dt * gravity(1);

    // #pragma omp parallel for collapse(2) num_threads(n_threads)
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            T mi = grid.mass(i,j);
            if (mi > 1e-25){

                T velocity_increment_x = -dt_particle_volume * grid_force_x(i,j) / mi + dt_gravity_x;
                T velocity_increment_y = -dt_particle_volume * grid_force_y(i,j) / mi + dt_gravity_y;

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
