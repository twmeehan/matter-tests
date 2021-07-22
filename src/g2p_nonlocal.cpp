#include "simulation.hpp"
#include <omp.h>


void Simulation::G2P_nonlocal(){

    std::fill( particles.delta_gamma_nonloc.begin(), particles.delta_gamma_nonloc.end(), 0.0 );
    std::vector<T> divisor; divisor.resize(Np); std::fill( divisor.begin(), divisor.end(), 0.0 );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<T> particles_delta_gamma_local(Np);
        std::vector<T> particles_divisor_local(Np);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            T delta_gamma_p = 0;
            T divisor_p = 0;
            unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - (nonlocal_support-1); // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - (nonlocal_support-1);
            unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - (nonlocal_support-1);

            for(int i = i_base; i < i_base+(2*nonlocal_support); i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+(2*nonlocal_support); j++){
                    T yi = grid.y[j];
                    for(int k = k_base; k < k_base+(2*nonlocal_support); k++){
                        T zi = grid.z[j];
                        T dist_norm_sq = (xp(0)-xi)*(xp(0)-xi) + (xp(1)-yi)*(xp(1)-yi) + (xp(2)-zi)*(xp(2)-zi);
                        if (dist_norm_sq <= nonlocal_l_sq){
                            T kernel = std::exp(-4.83597586204940892215090053992 * dist_norm_sq / nonlocal_l_sq);
                            T mi = grid.mass[ind(i,j,k)];
                            delta_gamma_p += grid.delta_gamma[ind(i,j,k)] * kernel * mi;
                            divisor_p     += kernel * mi;
                        } // end if
                    } // end loop k
                } // end loop j
            } // end loop i

            particles_delta_gamma_local[p] = delta_gamma_p;
            particles_divisor_local[p] = divisor_p;
        } // end loop p

        #pragma omp critical
        {
            for(int p = 0; p < Np; p++){
                particles.delta_gamma_nonloc[p] += particles_delta_gamma_local[p];
                divisor[p]                      += particles_divisor_local[p];
            }
        } // end omp critical

    } // end omp paralell

    for (int p = 0; p<Np; p++){
        T divisor_p = divisor[p];
        if (divisor_p > 0)
            particles.delta_gamma_nonloc[p] /= divisor_p;
        else
            particles.delta_gamma_nonloc[p] = 0;
        //grid.v[l] = (mi > 0) ? grid.v[l]/mi : TV::Zero(); // condition ? result_if_true : result_if_false
    }

} // end G2P_Optimized_Parallel
