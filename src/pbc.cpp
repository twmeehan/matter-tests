#include "simulation.hpp"

void Simulation::PBC(unsigned int safety_factor){

    for(int j = 0; j < Ny; j++){

        grid.v[ind(0,j)] = grid.v[ind(Nx-2-safety_factor-1,j)];
        grid.v[ind(1,j)] = grid.v[ind(Nx-1-safety_factor-1,j)];

        grid.v[ind(Nx-3,j)] = grid.v[ind(safety_factor,  j)];
        grid.v[ind(Nx-2,j)] = grid.v[ind(safety_factor+1,j)];
        grid.v[ind(Nx-1,j)] = grid.v[ind(safety_factor+2,j)];



        grid.flip[ind(0,j)] = grid.flip[ind(Nx-2-safety_factor-1,j)];
        grid.flip[ind(1,j)] = grid.flip[ind(Nx-1-safety_factor-1,j)];

        grid.flip[ind(Nx-3,j)] = grid.flip[ind(safety_factor,  j)];
        grid.flip[ind(Nx-2,j)] = grid.flip[ind(safety_factor+1,j)];
        grid.flip[ind(Nx-1,j)] = grid.flip[ind(safety_factor+2,j)];



        grid.mass[ind(0,j)] = grid.mass[ind(Nx-2-safety_factor-1,j)];
        grid.mass[ind(1,j)] = grid.mass[ind(Nx-1-safety_factor-1,j)];

        grid.mass[ind(Nx-3,j)] = grid.mass[ind(safety_factor,  j)];
        grid.mass[ind(Nx-2,j)] = grid.mass[ind(safety_factor+1,j)];
        grid.mass[ind(Nx-1,j)] = grid.mass[ind(safety_factor+2,j)];

    }


} // end pbc()
