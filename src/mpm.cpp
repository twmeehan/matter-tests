#include "simulation.hpp"
#include "tools.hpp"
#include "sampling_particles.hpp"

#include "vdb_obj.hpp" // Comment if not using openvdb

int main(){
    openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.directory = "/media/blatny/LaCie/larsie/vdb_test/";
    sim.sim_name = "lecture_dp_fric075";
    sim.end_frame = 200;
    sim.fps = 10;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 39;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1000); // feeder, incl slope and gran collapse

    // sim.Lx = 10;
    // sim.Ly = 10;
    // SampleParticles(sim.Lx, sim.Ly,    0.1, 6, 0, sim);
    // for(int p = 0; p < sim.Np; p++){
    //     sim.particles.x[p](0) += 100;
    //     sim.particles.x[p](1) += 48;
    // }

    sim.Lx = 0.3;
    sim.Ly = 0.15;
    T k_rad = 0.002;
    SampleParticles(sim.Lx, sim.Ly,     k_rad, 6, 5, sim);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    auto new_part_x = sim.particles.x;
    TV tmp_part(0.0, 0.0);
    new_part_x.push_back(tmp_part);
    sim.Np += 1;
    sim.particles = Particles(sim.Np);
    sim.particles.x = new_part_x;

    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;

    std::string name;
    // name = "Terrain"; VdbObj     terrain = VdbObj("/home/blatny/Documents/LevelSets/slope40degree_minus50.vdb", STICKY, 0.5, name); sim.objects.push_back(&terrain);
    name = "Ground"; PlateObj ground = PlateObj(0,  1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,     1,0);  sim.plates.push_back(ground);
    // name = "Ground";  AnalyticObj ground = AnalyticObj(SEPARATE, 0.5, name,   0, 0.43);   sim.objects_general.push_back(&ground);
    // name = "Gate";    AnalyticObj gate   = AnalyticObj(SEPARATE, 0,   name, 100, 0.016);  sim.objects_general.push_back(&gate);

    // TV test_point(0,1);
    // bool inside = terrain.inside(test_point);
    // debug("Is inside: ", inside);
    // TV normal = terrain.normal(test_point);
    // debug("Normal vector: \n", normal);

    // Elastoplasticity
    sim.elastic_model = StvkWithHencky;
    // sim.plastic_model = NoPlasticity;
    // sim.plastic_model = VonMises;
    sim.plastic_model = DruckerPrager;
    // sim.plastic_model = PerzynaVM;
    // sim.plastic_model = PerzynaDP;
    // sim.plastic_model = PerzynaMuIDP;
    // sim.plastic_model = MCC;
    // sim.plastic_model = MCCHard;
    // sim.plastic_model = MCCHardExp;
    // sim.plastic_model = PerzynaMCC;
    // sim.plastic_model = PerzynaMCCHard;
    // sim.plastic_model = PerzynaMuIMCC;
    // sim.plastic_model = PerzynaSinterMCC;

    sim.dp_slope = 1.5;
    sim.dp_cohesion = 0;

    sim.simulate();

	return 0;
}
