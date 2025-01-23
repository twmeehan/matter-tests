// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include <gtest/gtest.h>

#include "../src/simulation.hpp"
#include "../src/tools.hpp"

#include "../src/objects/object_vdb.hpp"

TEST(BoundaryTest, LoadFile) {
    std::string file_path = std::string(INCLUDE_DIR) + "/curve.vdb";
    std::ifstream file(file_path);
    ASSERT_TRUE(file.is_open()) << "Failed to open file at: " << file_path;
}

TEST(BoundaryTest, VDB) {

    openvdb::initialize();

    Simulation sim;
    sim.initialize(false);

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.gravity[1] = -9.81;

    sim.dx = 0.001;
    sim.particle_volume = std::pow(sim.dx, sim.dim);
    sim.particle_mass = sim.rho * sim.particle_volume;
    sim.Np = 1;
    sim.particles = Particles(sim.Np);

    sim.particles.x[0](0) = -1;
    sim.particles.x[0](1) = 1;

    std::string file_path = std::string(INCLUDE_DIR) + "/curve.vdb";
    std::string name = "Curve"; ObjectVdb curve  = ObjectVdb(file_path, SEPARATE, 0.0, name); sim.objects.push_back(&curve);


    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}
