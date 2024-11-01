#include "simulation.hpp"
#include "tools.hpp"
#include "sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_chute.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
#include "objects/object_vdb.hpp"
#include "sampling_particles_vdb.hpp"


int main(){
    openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.directory = "output/";
    sim.sim_name = "collapse";

    sim.end_frame = 20;     // Last frame to simulate
    sim.fps = 10;           // frames per second (float)
    sim.n_threads = 8;      // number of threads in parallel (integer)
    sim.save_grid = true;   // [default: false] save grid data
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].
    //    Positive numbers  [0,1]: PIC/FLIP where 1 equals pure FLIP
    //    Negative numbers [-1,0): APIC/AFLIP where -1 equals pure AFLIP

    sim.pbc = false; // [default: false] if true, use periodic boundary conditions, see pbc.cpp

    // INITILIZE ELASTICITY AND ELASTIC PARAMETERS
    sim.elastic_model = StvkWithHencky;
    sim.initialize(/* Young's (Pa) */ 1e6, /* Poisson's (-) */ 0.3, /* Density (kg/m3) */ 1000);
    sim.use_von_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)

    ////// GRAVITY ANGLE [default: gravity is 0]
    T theta_deg = 0; // angle in degrees of gravity vector, 0 means in negative y-direction
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);
    sim.gravity_special = false; // [default: false] if true, you can create a special gravity function in updateDt.cpp


    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad = 0.01;
    #ifdef THREEDIM
        sim.Lz = 5;
        SampleParticles(sim, sim.Lx, sim.Ly, sim.Lz, k_rad);
    #else
        SampleParticles(sim, sim.Lx, sim.Ly, k_rad);
    #endif

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

    ///// SET MAX TIME STEP
    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;

    ////// OBJECTS AND TERRAINS
    T friction = 0.2; // used if SEPARATE or SLIP
    std::string name;
    name = "Ground";  ObjectPlate ground = ObjectPlate(0,  1e10, -1e10, bottom, STICKY, friction, name);  sim.plates.push_back(ground);

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
    // name = "Bump";    ObjectBump bump    = ObjectBump(SEPARATE, friction, name);  sim.objects.push_back(&bump);
    // name = "Gate";    ObjectGate gate    = ObjectGate(SEPARATE, friction, name);  sim.objects.push_back(&gate);

    /////// Here is an example how to use ObjectVdb:
    // name = "Terrain"; ObjectVdb terrain  = ObjectVdb("../levelsets/vdb_file_name.vdb", STICKY, friction, name); sim.objects.push_back(&terrain);

    ////// OPTIONAL: TEST YOUR TERRAIN
    // TV test_point(0,1);
    // bool inside = terrain.inside(test_point);
    // debug("Is inside: ", inside);
    // TV normal = terrain.normal(test_point);
    // debug("Normal vector: \n", normal);

    ////// PLASTICITY MODEL
    ////// Choose ONE of the below models
    // sim.plastic_model = NoPlasticity;   // This must be chosen if pure elastic simulation
    // sim.plastic_model = VonMises;       // Von Mises
    // sim.plastic_model = DruckerPrager;  // Drucker-Prager
    // sim.plastic_model = DPSoft;         // Drucker-Prager with softening
    // sim.plastic_model = PerzynaVM;      // Perzyna model with von Mises yield
    sim.plastic_model = PerzynaDP;         // Perzyna model with Drucker yield
    // sim.plastic_model = PerzynaMuIDP;   // Drucker-Prager mu(I) rheology
    // sim.plastic_model = PerzynaMuIMCC;  // Modified cam clay mu(I) rheology
    // sim.plastic_model = MCC;            // Modified cam clay with explicit hardening
    // sim.plastic_model = MCCHardExp;     // Modified cam clay with implicit hardening
    // sim.plastic_model = PerzynaMCC;     // Perzyna model with MCC yield

    sim.use_pradhana = true; // [default: true] Use true to supress unwanted volume expansion in Drucker-Prager models

    sim.use_von_mises_q = false; // [default: false] if true, q is defined as q = sqrt(3/2 * s:s), otherwise q = sqrt(1/2 * s:s)
                                 // Note if using plastic model DPSoft then this will always be true!

    ////// PLASTIC PARAMETERS
    ////// Drucker-Prager models:
    sim.dp_slope = std::tan(30*M_PI/180.0); // [default: 1] Internal friction
    sim.dp_cohesion = 0; // [default: 0] Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case

    ////// Perzyna models:
    sim.perzyna_exp = 1; // [default: 1] Exponent in Perzyna models
    sim.perzyna_visc = 0; // [default: 0] Viscous time parameter is Perzyna models


    // Von Mises models:
    sim.yield_stress_orig = 100; // [default: 100]
    sim.yield_stress_min = 100; // [default: 100]
    sim.vm_ptensile = -1e20; // [default: -1e20]
    sim.xi = 0; // [default: 0]

     // Modified Cam-Clay models:
    sim.M = 1; // [default: 1] Slope of critical state line
    sim.beta = 0; // [default: 0] Cohesion (ratio of tensile strength to compressive strength)
    sim.p0 = 1000; // [default: 1000]

    // Mu(I) rheology models:
    sim.rho_s = sim.rho / 0.63; // Intrinsic grain density
    sim.grain_diameter = 1e-3; // Intrinsic grain diameter
    sim.in_numb_ref = 0.279;   // I_0
    sim.mu_1 = std::tan(20.9*M_PI/180.0); // mu_1
    sim.mu_2 = std::tan(32.8*M_PI/180.0); // mu_2

    sim.simulate();

	return 0;
}
