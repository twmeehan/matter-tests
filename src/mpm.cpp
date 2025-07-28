// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_box.hpp"
#include "objects/object_plate.hpp"
#include <random>
#include <cmath>
#include "./data_structures.hpp"




////////////////////////////////////////////////////////
// Setup environment for simulation
////////////////////////////////////////////////////////

void setup_environment_with_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::left,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::right,  BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::back,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::front,  BC::NoSlip));

    TV center(0, 1, 0);
    TV half_extents(0.4, 0.4, 0.4);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_without_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));

    TV center(0, 1, 0);
    TV half_extents(0.4, 0.4, 0.4);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_for_disney(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));

    TV center(0.5, 0.4, 0.5);
    TV half_extents(0.1, 0.1, 0.5); // half of (0.2, 0.2, 1.0)

    TM R = TM::Identity();         // Rotation matrix
    T mytheta = M_PI / 4;          // 45 degrees
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);
    // Z-axis rotation → top-down spin

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_only_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::left,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::right,  BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::back,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::front,  BC::NoSlip));
}

////////////////////////////////////////////////////////
// Setup shared simulation factors
////////////////////////////////////////////////////////

void setup_simulation_with_particles(int particle_count, Simulation& sim) {
    T theta_deg = 0;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad;
    switch (particle_count) {
        case 20000: k_rad = 0.031; break;
        case 10000: k_rad = 0.039; break;
        case 5000:  k_rad = 0.050; break;
        case 2000:  k_rad = 0.068; break;
        case 3000:  k_rad = 0.059; break;
        default:    k_rad = 0.05;  break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif
    T ppc = 8;
    sim.dx = std::cbrt(ppc / T(particle_count) * sim.Lx*sim.Ly*sim.Lz);
    sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc; // = Lx*Ly*Lz / T(square_samples.size())
    sim.particle_mass = sim.rho * sim.particle_volume;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    
    sim.Np = particle_count;
    sim.particles = Particles(sim.Np);
    for(int p = 0; p < sim.Np; p++){
        double x = dist(gen)*0.4;
        double y = dist(gen)*0.2;
        double z = dist(gen)*0.2;
        sim.particles.x[p](0)=x;
        sim.particles.x[p](1)=y;
        sim.particles.x[p](2)=z;
    }
    sim.grid_reference_point = TV::Zero();
}

////////////////////////////////////////////////////////
// Simulations
////////////////////////////////////////////////////////

int sand_collision(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "sand_collision_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_with_bounding_box(sim);

    sim.plastic_model = PlasticModel::DPVisc;
    sim.use_pradhana = true;
    sim.use_mises_q = false;
    sim.M = std::tan(30*M_PI/180.0);
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

int bouncy_cube(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "bouncy_cube_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = 0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e5;
    sim.nu = 0.45;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_with_bounding_box(sim);

    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.hardening_law = HardeningLaw::NoHard;
    sim.use_pradhana = false;
    sim.use_mises_q = false;
    sim.M = 0;
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 3000;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

int no_plasticity(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "no_plasticity_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_without_bounding_box(sim);

    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.use_pradhana = true;
    sim.use_mises_q = false;
    sim.M = std::tan(30*M_PI/180.0);
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

int snow(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "snow_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  += 0.3;
        sim.particles.x[p][1]  += 0.7;
        sim.particles.x[p][2]  += 0.4;
    }

    setup_environment_for_disney(sim);

    sim.plastic_model = PlasticModel::DPVisc;
    sim.use_pradhana = true;
    sim.use_mises_q = false;
    sim.M = std::tan(30*M_PI/180.0);
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////

void record_to_csv(const std::string& test_name, int particle_count, int duration_ms) {
    std::ofstream out("results.csv", std::ios::app);
    if (!out) {
        std::cerr << "Failed to open results.csv\n";
        return;
    }
    out << test_name << "," << particle_count << "," << duration_ms << "\n";
}

int main() {
    int duration;
    duration = snow(100);
    record_to_csv("snow", 100, duration);

    return 0;
}
