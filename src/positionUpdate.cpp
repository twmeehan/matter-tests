#include "simulation.hpp"


void Simulation::positionUpdate(){

    #pragma omp parallel for num_threads(n_threads)
    for(int p=0; p<Np; p++){

        //// Position is updated according to PIC velocities
        particles.x[p] = particles.x[p] + dt * particles.pic[p];

        //// Velicity is updated
        if (flip_ratio < -1){ // APIC
            particles.v[p] = particles.pic[p];
        } else if (flip_ratio < 0){ // AFLIP
            particles.v[p] = (-flip_ratio) * ( particles.v[p] + particles.flip[p] ) + (1 - (-flip_ratio)) * particles.pic[p];
        } else{ // PIC-FLIP
            particles.v[p] =   flip_ratio  * ( particles.v[p] + particles.flip[p] ) + (1 -   flip_ratio)  * particles.pic[p];
        }


        if (pbc){
            if (particles.x[p](0) > Lx){
                particles.x[p](0) = particles.x[p](0) - Lx;
            }
            else if (particles.x[p](0) < 0){
                particles.x[p](0) = Lx + particles.x[p](0);
            }
        }

        if (pbc_special){
            if (particles.x[p](0) > 0.8){ // 0.8 for nico bump
                particles.x[p](0)  = -0.1;
                particles.x[p](1) +=  0.27;
            }
            // if (particles.x[p](0) > 1.8){
            //     particles.x[p](0)  = -0.18;
            //     particles.x[p](1) +=  0.5;

               // particles.v[p](0) *= 0.1;
            // }
        }

    } // end loop over particles
}
