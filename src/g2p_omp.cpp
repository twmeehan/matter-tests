#include "simulation.hpp"
#include <omp.h>


void Simulation::G2P_Optimized_Parallel(){

    std::fill( particles.pic.begin(), particles.pic.end(), TV::Zero() );
    std::fill( particles.flip.begin(), particles.flip.end(), TV::Zero() );
    std::fill( particles.reg_laplacian.begin(), particles.reg_laplacian.end(), 0.0 );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> particles_pic_local(Np); std::fill( particles_pic_local.begin(), particles_pic_local.end(), TV::Zero() );
        std::vector<TV> particles_flip_local(Np); std::fill( particles_flip_local.begin(), particles_flip_local.end(), TV::Zero() );
        std::vector<T> particles_reg_local(Np);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            TV vp = TV::Zero();
            TV flipp = TV::Zero();
            T reg_laplacian_p = 0;
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
                        vp               += grid.v[ind(i,j,k)]    * weight;
                        flipp            += grid.flip[ind(i,j,k)] * weight;
                        reg_laplacian_p  += grid.reg_laplacian[ind(i,j,k)] * weight;
                    } // end loop k
                } // end loop j
            } // end loop i
            particles_pic_local[p]  = vp;
            particles_flip_local[p] = flipp;
            particles_reg_local[p]  = reg_laplacian_p;
        } // end loop p

        #pragma omp critical
        {
            for(int p = 0; p < Np; p++){
                particles.pic[p]            += particles_pic_local[p];
                particles.flip[p]           += particles_flip_local[p];
                particles.reg_laplacian[p] += particles_reg_local[p];
            }
        } // end omp critical

    } // end omp paralell

} // end G2P_Optimized_Parallel
