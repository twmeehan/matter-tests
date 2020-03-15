#include "simulation.hpp"
#include <omp.h>

void Simulation::P2G_Optimized_Parallel(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    grid.mass.setZero(Nx, Ny);
    grid.vx.setZero(Nx, Ny);
    grid.vy.setZero(Nx, Ny);

    #pragma omp parallel num_threads(n_threads)
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
