#include "simulation.hpp"

// Remember:
// P = dPsidF              (first Piola-Kirchhoff stress tensor)
// tau = P * F.transpose() (Kirchhoff stress tensor)
/*
void Simulation::explicitEulerUpdate_Baseline(){
    TV grad_wip;
    TM Fe, dPsidF;

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    // TV external_gravity;
    ////////////////////////////////////////////////////////

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if (grid.mass(i,j) > 1e-25){
                T xi = grid.x(i);
                T yi = grid.y(j);
                TV grid_force = TV::Zero();
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
                            Eigen::JacobiSVD<TM> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
                            TA sigma = svd.singularValues().array(); // abs() for inverse also??
                            TM logSigma = sigma.abs().log().matrix().asDiagonal();
                            TM invSigma = sigma.inverse().matrix().asDiagonal();
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

                TV velocity_increment = -dt * particle_volume * grid_force / grid.mass(i,j) + dt * gravity;

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
*/

void Simulation::explicitEulerUpdate_Optimized(){
    TM Fe, dPsidF, tau;

    std::vector<TV> grid_force(Nx*Ny*Nz, TV::Zero());

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

        TV xp = particles.x[p];
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1;
        unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - 1;

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x[i];
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y[j];
                for(int k = k_base; k < k_base+4; k++){
                    T zi = grid.z[k];
                    if (grid.mass[ind(i,j,k)] > 0){
                        grid_force[ind(i,j,k)] += tau * grad_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                    } // end if non-zero grid mass
                } // end for k
             } // end for j
         } // end for i
    } // end for particles

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    ////////////////////////////////////////////////////////

    T dt_particle_volume = dt * particle_volume;
    TV dt_gravity = dt * gravity;

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            for(int k = 0; k < Nz; k++){
                T mi = grid.mass[ind(i,j,k)];
                if (mi > 0){

                    TV velocity_increment = -dt_particle_volume * grid_force[ind(i,j,k)] / mi + dt_gravity;

                    //////////// if external grid gravity: //////////////////
                    // T external_gravity = external_gravity_pair.first(i,j);
                    // T external_gravity = external_gravity_pair.second(i,j);
                    // velocity_increment_x += dt * external_gravity(0);
                    // velocity_increment_y += dt * external_gravity(1);
                    ////////////////////////////////////////////////////////

                    TV old_vi = grid.v[ind(i,j,k)];
                    TV new_vi = old_vi + velocity_increment;
                    boundaryCollision(grid.x[i], grid.y[j], grid.z[k], new_vi);

                    // Currently not working:
                    // boundaryCorrection(new_xi, new_yi, new_vxi, new_vyi);

                    // Only if impose velocity on certain grid nodes:
                    // overwriteGridVelocity(grid.x[i], grid.y[j], new_vi);

                    grid.v[ind(i,j,k)] = new_vi;
                    grid.flip[ind(i,j,k)] = new_vi - old_vi;

                } // end if non-zero grid mass
            } // end for k
        } // end for j
    } // end for i

} // end explicitEulerUpdate_Optimized
