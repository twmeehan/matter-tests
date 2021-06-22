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

    grid.v.resize(Nx*Ny*Nz); std::fill( grid.v.begin(), grid.v.end(), TV::Zero() );
    grid.mass.resize(Nx*Ny*Nz); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );
    grid.reg_laplacian.resize(Nx*Ny*Nz); std::fill( grid.reg_laplacian.begin(), grid.reg_laplacian.end(), 0.0 );

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
                        grid.mass[ind(i,j,k)]           += weight;
                        grid.v[ind(i,j,k)]              += particles.v[p] * weight;
                        grid.reg_laplacian[ind(i,j,k)] += particles.reg_variable[p] * laplace_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx, one_over_dx_square);
                    }
                } // end for k
            } // end for j
        } // end for i
    } // end for p

    // Apply BC to regularization here!

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

} // end P2G_Optimized
