// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include "tools.hpp"
#include "sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_chute.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
// #include "objects/object_vdb.hpp"
// #include "sampling_particles_vdb.hpp"


int main(){
    // openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "collapse");
    sim.save_grid = true;   // save grid data

    sim.end_frame = 20;     // last frame to simulate
    sim.fps = 10;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].
    sim.reduce_verbose = true; // reduce the screen output

    // INITILIZE ELASTICITY
    sim.elastic_model = Hencky;
    sim.E = 1e6;     // Young's modulus (Pa)
    sim.nu = 0.3;   // Poisson's ratio (-)
    sim.rho = 1000; // Density (kg/m3)

    ////// GRAVITY ANGLE [default: gravity is 0]
    T theta_deg = 0; // angle in degrees of gravity vector
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad = 0.01;
    #ifdef THREEDIM
        sim.Lz = 5;
    #endif
    SampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    ////// OPTIONAL: ADD INDIVIDUAL PARTICLES
    auto new_particle_x = sim.particles.x;
    #ifdef THREEDIM
        TV tmp_particle(0.0, 0.0, 0.0);
    #else
        TV tmp_particle(0.0, 0.0);
    #endif
    new_particle_x.push_back(tmp_particle);
    sim.Np += 1;
    sim.particles = Particles(sim.Np);
    sim.particles.x = new_particle_x;

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
    T friction = 0.2; // used if SEPARATE or SLIP
    std::string name;
    name = "Ground";  ObjectPlate ground = ObjectPlate(0,  1e10, -1e10, bottom, STICKY, friction, name);  sim.plates.push_back(ground);

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
    // name = "Bump";    ObjectBump bump    = ObjectBump(SEPARATE, friction, name);  sim.objects.push_back(&bump);
    // name = "Gate";    ObjectGate gate    = ObjectGate(SEPARATE, friction, name);  sim.objects.push_back(&gate);

    /////// Here is an example how to use ObjectVdb:
    // name = "Terrain"; ObjectVdb terrain  = ObjectVdb("../levelsets/vdb_file_name.vdb", STICKY, friction, name); sim.objects.push_back(&terrain);

    ////// PLASTICITY
    sim.plastic_model = DPVisc; // Perzyna model with Drucker_Prager yield surface

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.use_von_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    sim.dp_slope = std::tan(30*M_PI/180.0); // Internal friction
    sim.dp_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 0; // Viscous time parameter is Perzyna models

    sim.simulate();

	return 0;
}
