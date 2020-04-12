#include "simulation.hpp"
#include <omp.h>

void Simulation::P2G_Optimized_Parallel(){

    grid.v.resize(Nx*Ny*Nz); std::fill( grid.v.begin(), grid.v.end(), TV::Zero() );
    grid.mass.resize(Nx*Ny*Nz); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );
    grid.regularization.resize(Nx*Ny*Nz); std::fill( grid.regularization.begin(), grid.regularization.end(), 0.0 );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> grid_v_local(Nx*Ny*Nz, TV::Zero() );
        std::vector<T> grid_mass_local(Nx*Ny*Nz);
        std::vector<T> grid_reg_local(Nx*Ny*Nz);

        #pragma omp for
        for(int p = 0; p < Np; p++){
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
                        T weight = wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);

                        if (weight > 1e-25){
                            grid_mass_local[ind(i,j,k)] += weight;
                            grid_v_local[ind(i,j,k)]    += particles.v[p] * weight;
                            grid_reg_local[ind(i,j,k)]  += particles.eps_pl_dev[p] * laplace_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx, one_over_dx_square);
                        }
                    } // end for k
                } // end for j
            } // end for i
        } // end for p


        #pragma omp critical
        {
            for (int l = 0; l<Nx*Ny*Nz; l++){
                grid.mass[l]           += grid_mass_local[l];
                grid.v[l]              += grid_v_local[l];
                grid.regularization[l] += grid_reg_local[l];
            } // end for l
        } // end omp critical

    } // end omp parallel

    ///////////////////////////////////////////////////////////
    // At this point in time grid.mass is equal to m_i / m_p //
    ///////////////////////////////////////////////////////////

    for (int l = 0; l<Nx*Ny*Nz; l++){
        T mi = grid.mass[l];
        if (mi > 0)
            grid.v[l] /= mi;
        else
            grid.v[l].setZero();
        //grid.v[l] = (mi > 0) ? grid.v[l]/mi : TV::Zero(); // condition ? result_if_true : result_if_false
    }

    for(auto&& m: grid.mass){
        m *= particle_mass;
    }

} // end P2G_Optimized_Parallel
