#include "simulation.hpp"
#include <omp.h>


void Simulation::explicitEulerUpdate_Optimized_Parallel(){

    std::vector<TV> grid_force(Nx*Ny*Nz, TV::Zero());

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> grid_force_local(Nx*Ny*Nz, TV::Zero());

        #pragma omp for
        for(int p = 0; p < Np; p++){

            TM Fe = particles.F[p];

            TM dPsidF;
            if (elastic_model == NeoHookean){
                dPsidF = NeoHookeanPiola(Fe);
            }
            else if (elastic_model == StvkWithHencky){ // St Venant Kirchhoff with Hencky strain
                dPsidF = StvkWithHenckyPiola(Fe);
            }
            else{
                debug("You specified an unvalid ELASTIC model!");
            }

            TM tau = dPsidF * Fe.transpose();
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
                        if ( grid.mass[ind(i,j,k)] > 0){
                            grid_force_local[ind(i,j,k)] += tau * grad_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                        } // end if non-zero grid mass
                    } // end for k
                 } // end for j
             } // end for i
        } // end for particles

        #pragma omp critical
        {
            for (int l = 0; l<Nx*Ny*Nz; l++){
                grid_force[l] += grid_force_local[l];
            } // end for l
        } // end omp critical

    } // end omp parallel

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    ////////////////////////////////////////////////////////

    T dt_particle_volume = dt * particle_volume;
    TV dt_gravity = dt * gravity;

    #pragma omp parallel for collapse(3) num_threads(n_threads)
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

} // end explicitEulerUpdate_Optimized_Parallel
