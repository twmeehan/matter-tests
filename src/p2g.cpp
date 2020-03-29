#include "simulation.hpp"

/*
void Simulation::P2G_Baseline(){
    // This (nested) loop over i (and j) can easily be paralellized
    // Need to create thread-local grid.mass, grid.vx, grid.vy
    for(int i=0; i<Nx; i++){
        T xi = grid.x(i);
        for(int j=0; j<Ny; j++){
            T yi = grid.y(j);
            T vxi = 0;
            T vyi = 0;
            T mass = 0;
            for(int p=0; p<Np; p++){
                T xp = particles.x(p);
                T yp = particles.y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(yp-yi) < 1.5*dx){
                    T weight = wip(xp, yp, xi, yi, one_over_dx);
                    mass += weight;
                    vxi  += particles.vx(p) * weight;
                    vyi  += particles.vy(p) * weight;
                }
            } // end for particles
            grid.mass(i,j) = mass * particle_mass;
            if (mass < 1e-25){
                grid.vx(i,j)   = 0.0;
                grid.vy(i,j)   = 0.0;
            } else {
                grid.vx(i,j)   = vxi / mass;
                grid.vy(i,j)   = vyi / mass;
            }
        } // end for j
    } // end for i
} // end P2G_Baseline
*/

void Simulation::P2G_Optimized(){

    grid.vx.setZero(Nx, Ny);
    grid.vy.setZero(Nx, Ny);
    grid.regularization.setZero(Nx, Ny);

    for(int p = 0; p < Np; p++){
        TV2 xp = particles.x[p];
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);
                T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);

                if (weight > 1e-25){
                    grid.mass(i,j) += weight;
                    grid.vx(i,j)   += particles.v[p](0) * weight;
                    grid.vy(i,j)   += particles.v[p](1) * weight;
                    grid.regularization(i,j) += particles.eps_pl_dev[p] * laplace_wip(xp(0), xp(1), xi, yi, one_over_dx, one_over_dx_square);
                }

            } // end for j
        } // end for i
    } // end for p

    // Apply BC to regularization here!

    ///////////////////////////////////////////////////////////
    // At this point in time grid.mass is equal to m_i / m_p //
    ///////////////////////////////////////////////////////////
    grid.vx = (grid.mass.array() > 0).select( (grid.vx.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.vy = (grid.mass.array() > 0).select( (grid.vy.array() / grid.mass.array()).matrix(), TMX::Zero(Nx, Ny) );
    grid.mass *= particle_mass;

} // end P2G_Optimized
