#include "simulation.hpp"
#include "tools.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_chute.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
#include "objects/object_vdb.hpp"
#include "sampling/sampling_particles_from_vdb.hpp"


int main(){
    openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.directory = "output/";
    sim.sim_name = "my_simulation_name";

    sim.end_frame = 20;     // Last frame to simulate
    sim.fps = 10;           // frames per second (float)
    sim.n_threads = 8;      // number of threads in parallel (integer)
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].
    //    Positive numbers  [0,1]: PIC/FLIP where 1 equals pure FLIP
    //    Negative numbers [-1,0): APIC/AFLIP where -1 equals pure AFLIP

    // INITILIZE ELASTICITY AND ELASTIC PARAMETERS
    sim.elastic_model = StvkWithHencky;
    sim.initialize(/* Young's (Pa) */ 1e6, /* Poisson's (-) */ 0.3, /* Density (kg/m3) */ 1000);

    ////// GRAVITY ANGLE
    T theta_deg = 0; // angle in degrees of gravity vector (0 means in negative y-direction)
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

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
    TV tmp_particle(0.0, 0.0);
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
    #ifdef THREEDIM
        name = "Ground";  ObjectPlate ground = ObjectPlate(0,  1e10, -1e10, bottom, STICKY, friction, name);  sim.plates.push_back(ground);
    #else
        name = "Ground";  ObjectPlate ground = ObjectPlate(0,  1e10, -1e10, bottom, STICKY, friction, name);  sim.plates.push_back(ground);
    #endif
    // name = "Bump";    ObjectBump bump    = ObjectBump(SEPARATE, friction, name);  sim.objects.push_back(&bump);
    // name = "Gate";    ObjectGate gate    = ObjectGate(SEPARATE, friction, name);  sim.objects.push_back(&gate);
    // name = "Terrain"; ObjectVdb terrain  = ObjectVdb("/absolute_path_to_directory/vdb_file_name.vdb", STICKY, friction, name); sim.objects.push_back(&terrain);

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

    sim.use_von_mises_q = false;
    sim.use_pradhana = true;

    ////// PLASTIC PARAMETERS
    sim.dp_slope = std::tan(30*M_PI/180.0);
    sim.dp_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

	return 0;
}
