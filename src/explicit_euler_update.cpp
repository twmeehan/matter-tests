#include "simulation.hpp"

// Remember:
// P = dPsidF              (first Piola-Kirchhoff stress tensor)
// tau = P * F.transpose() (Kirchhoff stress tensor)

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

                        grad_wip(0) = gradx_wip(xp, yp, xi, yi, one_over_dx);
                        grad_wip(1) = grady_wip(xp, yp, xi, yi, one_over_dx);

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

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    TMX grid_force_x = TMX::Zero(Nx, Ny);
    TMX grid_force_y = TMX::Zero(Nx, Ny);

    for(int p = 0; p < Np; p++){

        Fe = particles.F[p];

        if (elastic_model == NeoHookean){
            dPsidF = NeoHookeanPiola(Fe);
        }
        else if (elastic_model == StvkWithHencky){ // St Venant Kirchhoff with Hencky strain
            dPsidF = StvkWithHenckyPiola(Fe);
        }
        else{
            debug("You specified an unvalid ELASTIC model!");
        }

        tau = dPsidF * Fe.transpose();
        particles.tau[p] = tau;

        T xp = particles.x(p);
        T yp = particles.y(p);
        unsigned int i_base = std::floor((xp-x0)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);

                if (grid.mass(i,j) > 1e-25){
                    grad_wip(0) = gradx_wip(xp, yp, xi, yi, one_over_dx);
                    grad_wip(1) = grady_wip(xp, yp, xi, yi, one_over_dx);
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

                // Not working:
                // boundaryCorrection(new_xi, new_yi, new_vxi, new_vyi);

                // Only if impose velocity on certain grid nodes:
                // overwriteGridVelocity(new_xi, new_yi, new_vxi, new_vyi);

                grid.vx(i,j) = new_vxi;
                grid.vy(i,j) = new_vyi;

            } // end if non-zero grid mass

        } // end for j
    } // end for i

} // end explicitEulerUpdate_Optimized
