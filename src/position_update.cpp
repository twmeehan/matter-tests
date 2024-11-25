// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"


void Simulation::positionUpdate(){

    #pragma omp parallel for num_threads(n_threads)
    for(int p=0; p<Np; p++){

        //// Position is updated according to PIC velocities
        particles.x[p] = particles.x[p] + dt * particles.pic[p];

        if (musl == false){
            //// Velicity is updated
            if (flip_ratio < -1){ // APIC
                particles.v[p] = particles.pic[p];
            } else if (flip_ratio < 0){ // AFLIP
                particles.v[p] = (-flip_ratio) * ( particles.v[p] + particles.flip[p] ) + (1 - (-flip_ratio)) * particles.pic[p];
            } else{ // PIC-FLIP
                particles.v[p] =   flip_ratio  * ( particles.v[p] + particles.flip[p] ) + (1 -   flip_ratio)  * particles.pic[p];
            }
        }

        if (pbc){
            if (particles.x[p](0) > Lx){
                particles.x[p](0) = particles.x[p](0) - Lx;
            }
            else if (particles.x[p](0) < 0){
                particles.x[p](0) = Lx + particles.x[p](0);
            }
        }

        // To be hard-coded depending on application
        if (pbc_special){
            if (particles.x[p](0) > 0.8){
                particles.x[p](0)  = -0.1;
                particles.x[p](1) +=  0.27;
            }
        }

    } // end loop over particles
}
