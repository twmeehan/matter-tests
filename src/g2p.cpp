#include "simulation.hpp"

void Simulation::G2P_Baseline(){
    // This loop over p can be easily paralellized
    // Need to create thread-local versions of particles.vx and vy
    for(int p=0; p<Np; p++){
        T xp = particles.x(p);
        T yp = particles.y(p);
        T vxp = 0;
        T vyp = 0;
        for(int i=0; i<Nx; i++){
            T xi = grid.x(i);
            for(int j=0; j<Ny; j++){
                T yi = grid.y(j);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid.vx(i,j) * weight;
                    vyp += grid.vy(i,j) * weight;
                }
            }
        }
        particles.vx(p) = vxp;
        particles.vy(p) = vyp;
    }
} // end G2P_Baseline


void Simulation::G2P_Optimized(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    // This loop over p can be easily paralellized
    // Need to create thread-local versions of particles.vx and vy
    for(int p = 0; p < Np; p++){
        T xp = particles.x(p);
        T yp = particles.y(p);
        T vxp = 0;
        T vyp = 0;
        unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);
                T weight = wip(xp, yp, xi, yi, dx);
                vxp += grid.vx(i,j) * weight;
                vyp += grid.vy(i,j) * weight;
            } // end loop j
        } // end loop i
        particles.vx(p) = vxp;
        particles.vy(p) = vyp;
    } // end loop p
} // end G2P_Optimized

void Simulation::G2P_Optimized_Parallel(){

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    particles.vx.setZero(Np);
    particles.vy.setZero(Np);

    #pragma omp parallel num_threads(n_threads)
    {
        TVX particles_vx_local = TVX::Zero(Np);
        TVX particles_vy_local = TVX::Zero(Np);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            T xp = particles.x(p);
            T yp = particles.y(p);
            T vxp = 0;
            T vyp = 0;
            unsigned int i_base = std::floor((xp-x0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((yp-y0)/dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x(i);
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y(j);
                    T weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid.vx(i,j) * weight;
                    vyp += grid.vy(i,j) * weight;
                } // end loop j
            } // end loop i
            particles_vx_local(p) = vxp;
            particles_vy_local(p) = vyp;
        } // end loop p

        #pragma omp critical
        {
            for(int p = 0; p < Np; p++){
                particles.vx(p) += particles_vx_local(p);
                particles.vy(p) += particles_vy_local(p);
            }
        } // end omp critical

    } // end omp paralell

} // end G2P_Optimized_Parallel
