#include "simulation.hpp"
/*
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
                    T weight = wip(xp, yp, xi, yi, one_over_dx);
                    vxp += grid.vx(i,j) * weight;
                    vyp += grid.vy(i,j) * weight;
                }
            }
        }
        particles.vx(p) = vxp;
        particles.vy(p) = vyp;
    }
} // end G2P_Baseline
*/

void Simulation::G2P_Optimized(){

    // This loop over p can be easily paralellized
    // Need to create thread-local versions of particles.vx and vy
    for(int p = 0; p < Np; p++){
        TV2 xp = particles.x[p];
        TV2 vp = TV2::Zero();
        TV2 flipp = TV2::Zero();
        T regularization_p = 0;
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);
                T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                vp(0) += grid.vx(i,j) * weight;
                vp(1) += grid.vy(i,j) * weight;
                flipp(0) += grid.flipx(i,j) * weight;
                flipp(1) += grid.flipy(i,j) * weight;
                regularization_p += grid.regularization(i,j) * weight;
            } // end loop j
        } // end loop i
        particles.pic[p] = vp;
        particles.flip[p] = flipp;
        particles.regularization[p] = regularization_p;
    } // end loop p
} // end G2P_Optimized
