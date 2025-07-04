// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
<<<<<<< HEAD
=======
#include "objects/object_box.hpp"
>>>>>>> ed4816d (added my own tests)
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
// #include "objects/object_vdb.hpp"
// #include "sampling/sampling_particles_vdb.hpp"


<<<<<<< HEAD
int main(){
    // openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "collapse");

    sim.save_grid = true;
    sim.end_frame = 20;     // last frame to simulate
    sim.fps = 10;           // frames per second
=======
int sand_sim(int particle_count) {

    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    

    Simulation sim;
    std::string name = "sand_" + std::to_string(particle_count);

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ name);

    sim.save_grid = true;
    sim.end_frame = 60;     // last frame to simulate
    sim.fps = 30;           // frames per second
>>>>>>> ed4816d (added my own tests)
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
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
<<<<<<< HEAD
    T k_rad = 0.01;
    #ifdef THREEDIM
        sim.Lz = 5;
    #endif
=======
    T k_rad; //0.031 ---- 0.039 ---- 0.05 -- 0.068 -- 0.059
    switch (particle_count) {
        case 20000:
            k_rad = 0.031;
            break;
        case 10000:
            k_rad = 0.039;
            break;
        case 5000:
            k_rad = 0.050;
            break;
        case 2000:
            k_rad = 0.068;
            break;
        case 3000:
            k_rad = 0.059;
            break;
        default:
            k_rad = 0.05;
            break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif

>>>>>>> ed4816d (added my own tests)
    sampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
<<<<<<< HEAD
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
=======
        sim.particles.x[p](0) -= 0.5;
        sim.particles.x[p](1) += 3.0;
        sim.particles.x[p](2) -= 0.5;

>>>>>>> ed4816d (added my own tests)
    }
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
<<<<<<< HEAD
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip)); 

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
    // T friction = 0.2; 
    // sim.objects.push_back(std::make_unique<ObjectBump>(BC::SlipFree, friction));
    // sim.objects.push_back(std::make_unique<ObjectGate>(BC::SlipFree, friction));
=======
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

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
>>>>>>> ed4816d (added my own tests)

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
<<<<<<< HEAD
    sim.plastic_model = PlasticModel::DPVisc; // Perzyna model with Drucker_Prager yield surface
=======
    sim.plastic_model = PlasticModel::DPVisc;
>>>>>>> ed4816d (added my own tests)

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.use_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    sim.M = std::tan(30*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 0; // Viscous time parameter is Perzyna models

    sim.simulate();

<<<<<<< HEAD
	return 0;
=======
    auto end = high_resolution_clock::now(); // End timing
    auto duration = duration_cast<milliseconds>(end - start);
    return static_cast<int>(duration.count());
}

int bouncy_cube(int particle_count) {

    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    

    Simulation sim;
    std::string name = "test_" + std::to_string(particle_count);

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ name);

    sim.save_grid = true;
    sim.end_frame = 60;     // last frame to simulate
    sim.fps = 30;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = 0.95; // (A)PIC-(A)FLIP ratio in [-1,1].

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e5;     // Young's modulus (Pa)
    sim.nu = 0.45;   // Poisson's ratio (-)
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
    T k_rad; //0.031 ---- 0.039 ---- 0.05 -- 0.068 -- 0.059
    switch (particle_count) {
        case 20000:
            k_rad = 0.031;
            break;
        case 10000:
            k_rad = 0.039;
            break;
        case 5000:
            k_rad = 0.050;
            break;
        case 2000:
            k_rad = 0.068;
            break;
        case 3000:
            k_rad = 0.059;
            break;
        default:
            k_rad = 0.05;
            break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif

    sampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5;
        sim.particles.x[p](1) += 3.0;
        sim.particles.x[p](2) -= 0.5;

    }
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
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

    /////// Here are some examples how to use the objects derived from ObjectGeneral:

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
    sim.elastic_model     = ElasticModel::Hencky;   // Elasticity for large deformation
    sim.plastic_model     = PlasticModel::NoPlasticity;       // Perzyna viscoplastic model
    sim.hardening_law   = HardeningLaw::NoHard;       // No hardening
    sim.use_pradhana      = false;                      // Not needed for liquids
    sim.use_mises_q       = false;// [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    sim.M = 0; // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 3000; // Viscous time parameter is Perzyna models

    sim.simulate();

    auto end = high_resolution_clock::now(); // End timing
    auto duration = duration_cast<milliseconds>(end - start);
    return static_cast<int>(duration.count());
}

int test(int particle_count) {

    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    

    Simulation sim;
    std::string name = "test_" + std::to_string(particle_count);

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ name);

    sim.save_grid = true;
    sim.end_frame = 60;     // last frame to simulate
    sim.fps = 30;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = 0.95; // (A)PIC-(A)FLIP ratio in [-1,1].

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
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
    T k_rad; //0.031 ---- 0.039 ---- 0.05 -- 0.068 -- 0.059
    switch (particle_count) {
        case 20000:
            k_rad = 0.031;
            break;
        case 10000:
            k_rad = 0.039;
            break;
        case 5000:
            k_rad = 0.050;
            break;
        case 2000:
            k_rad = 0.068;
            break;
        case 3000:
            k_rad = 0.059;
            break;
        default:
            k_rad = 0.05;
            break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif

    sampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5;
        sim.particles.x[p](1) += 3.0;
        sim.particles.x[p](2) -= 0.5;

    }
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
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

    /////// Here are some examples how to use the objects derived from ObjectGeneral:

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
   // sim.elastic_model = ElasticModel::Hencky;

    sim.plastic_model = PlasticModel::DPVisc;

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.use_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    sim.M = std::tan(1200*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 0; // Viscous time parameter is Perzyna models

    sim.simulate();

    auto end = high_resolution_clock::now(); // End timing
    auto duration = duration_cast<milliseconds>(end - start);
    return static_cast<int>(duration.count());
}

int no_plasticity(int particle_count) {

    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    

    Simulation sim;
    std::string name = "no_plasticity_" + std::to_string(particle_count);

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ name);

    sim.save_grid = true;
    sim.end_frame = 60;     // last frame to simulate
    sim.fps = 30;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
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
    T k_rad; //0.031 ---- 0.039 ---- 0.05 -- 0.068 -- 0.059
    switch (particle_count) {
        case 20000:
            k_rad = 0.031;
            break;
        case 10000:
            k_rad = 0.039;
            break;
        case 5000:
            k_rad = 0.050;
            break;
        case 2000:
            k_rad = 0.068;
            break;
        case 3000:
            k_rad = 0.059;
            break;
        default:
            k_rad = 0.05;
            break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif

    sampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5;
        sim.particles.x[p](1) += 2.0;
        sim.particles.x[p](2) -= 0.5;

    }
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));

    TV center(0, 1, 0);
    TV half_extents(0.4, 0.4, 0.4);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));

    /////// Here are some examples how to use the objects derived from ObjectGeneral:

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.use_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    sim.M = std::tan(30*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 0; // Viscous time parameter is Perzyna models

    sim.simulate();

    auto end = high_resolution_clock::now(); // End timing
    auto duration = duration_cast<milliseconds>(end - start);
    return static_cast<int>(duration.count());
}

int main() {
    // openvdb::initialize(); // if needed

    int time_ms = test(5000);
    std::cout << "Sim took " << time_ms << " ms." << std::endl;

    return 0;
>>>>>>> ed4816d (added my own tests)
}
