// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

int main(){

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "cube_rotating");
    sim.save_grid = true;
    sim.reduce_verbose = true;

    sim.end_frame = 20;
    sim.fps = 1;

    sim.gravity = TV::Zero();

    sim.cfl = 0.5;
    sim.flip_ratio = -1;
    sim.n_threads = 8;

    T h_gate, l_gate;
    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    sampleParticles(sim, 0.01);

    T total_energy_init = 0;
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) -= 0.5*sim.Ly;

        T vx = -1.0*sim.particles.x[p](1) + 0.5;
        T vy =  1.0*sim.particles.x[p](0) + 0.5;
        sim.particles.v[p](0) = vx;
        sim.particles.v[p](1) = vy;

        total_energy_init += 0.5*(vx*vx + vy*vy); // per unit mass
    }

    // Elasticity
    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6;     // Young's modulus (Pa)
    sim.nu = 0.3;   // Poisson's ratio (-)
    sim.rho = 1550; // Density (kg/m3)

    sim.simulate();

    T total_energy_last = 0;
    for(int p = 0; p < sim.Np; p++){
        T vx = sim.particles.v[p](0);
        T vy = sim.particles.v[p](1);
        #ifdef THREEDIM
            T vz = sim.particles.v[p](2);
            total_energy_last += 0.5*(vx*vx + vy*vy + vz*vz); // per unit mass
        #else
            total_energy_last += 0.5*(vx*vx + vy*vy); // per unit mass
        #endif
    }

    debug("E_init = ", total_energy_init);
    debug("E_last = ", total_energy_last);

    T rel_diff = (total_energy_init - total_energy_last) / total_energy_init;
    debug("rel_diff = ", rel_diff);

	return 0;
}
