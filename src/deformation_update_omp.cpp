#include "simulation.hpp"
#include <omp.h>

// Deformation gradient is updated based on the NEW GRID VELOCITIES and the OLD PARTICLE POSITIONS
void Simulation::deformationUpdate_Parallel(){

    #ifdef WARNINGS
        debug("deformationUpdate_Parallel");
    #endif

    std::fill( particles.delta_gamma.begin(), particles.delta_gamma.end(), 0.0 );
    unsigned int plastic_count = 0;

    #pragma omp parallel for reduction(+:plastic_count) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        TM sum = TM::Zero();
        TV xp = particles.x[p];
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1;
    #ifdef THREEDIM
        unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - 1;
    #endif

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x[i];
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y[j];
    #ifdef THREEDIM
                for(int k = k_base; k < k_base+4; k++){
                    T zi = grid.z[k];
                    sum += grid.v[ind(i,j,k)] * grad_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx).transpose();
                } // end loop k
    #else
                sum += grid.v[ind(i,j)] * grad_wip(xp(0), xp(1), xi, yi, one_over_dx).transpose();
    #endif
            } // end loop i
        } // end loop j

        TM Fe_trial = particles.F[p];
        Fe_trial = Fe_trial + dt * sum * Fe_trial;
        particles.F[p] = Fe_trial;

        plasticity(p, plastic_count, Fe_trial);

    } // end loop over particles

    debug("               Proj: ", plastic_count, " / ", Np);

} // end deformationUpdate_Parallel
