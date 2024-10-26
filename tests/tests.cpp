#include <gtest/gtest.h>

#include "../src/simulation.hpp"
#include "../src/tools.hpp"
#include "../src/sampling_particles.hpp"




TEST(ElasticityTest, BulkModulus) {

    Simulation sim;
    sim.save_sim = false;
    sim.end_frame = 100;
    sim.fps = 1;
    sim.n_threads = 8;
    sim.cfl = 0.6;
    sim.flip_ratio = -0.95;

    sim.gravity = TV::Zero();

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 10000); // feeder, incl slope and gran collapse

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
    sim.Lz = 0.2;
        SampleParticles(sim, sim.Lx, sim.Ly, sim.Lz, 0.02, 8);
    #else
        SampleParticles(sim, sim.Lx, sim.Ly,         0.02, 4);
    #endif

    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;
    T vel = 0.001;

    T vmin_factor = 10;
    T load_factor = 1;

    std::string name;
    #ifdef THREEDIM
        name = "Ground";  ObjectPlate ground = ObjectPlate(0-0.5*sim.dx,       1e20, -1e20, bottom, STICKY, 0, name,   0,  vel, 0, vmin_factor, load_factor);  sim.plates.push_back(ground);
        name = "Compre";  ObjectPlate compre = ObjectPlate(sim.Ly+0.5*sim.dx,  1e20, -1e20, top,    STICKY, 0, name,   0, -vel, 0, vmin_factor, load_factor);  sim.plates.push_back(compre);
    #else
        name = "Ground";  ObjectPlate ground = ObjectPlate(0-0.5*sim.dx,       1e20, -1e20, bottom, STICKY, 0, name,   0,  vel,    vmin_factor, load_factor);  sim.plates.push_back(ground);
        name = "Compre";  ObjectPlate compre = ObjectPlate(sim.Ly+0.5*sim.dx,  1e20, -1e20, top,    STICKY, 0, name,   0, -vel,    vmin_factor, load_factor);  sim.plates.push_back(compre);
    #endif

    sim.elastic_model = StvkWithHencky;
    sim.plastic_model = NoPlasticity;

    sim.simulate();

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
    T Javg;
    sim.computeAvgData(volavg_cauchy, volavg_kirchh, Javg);

    T volavg_p = -1.0 * (volavg_kirchh(0,0) + volavg_kirchh(1,1)) / sim.dim;
    T volavg_epsv = std::log(Javg);

    T measured_K = volavg_p / (-volavg_epsv);

    T rel_diff = std::abs(measured_K - sim.K) / sim.K;

    ASSERT_NEAR(rel_diff, 0.0, 0.03);
}






TEST(GranularcollapseTest, DruckerPrager) {

    Simulation sim;
    sim.save_sim = false;
    sim.end_frame = 70;
    sim.fps = 50;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 10;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1000);

    sim.Lx = 0.20;
    sim.Ly = 0.15;
    T k_rad = 0.0015;
    #ifdef THREEDIM
        sim.Lz = 0.10;
        SampleParticles(sim, sim.Lx, sim.Ly, sim.Lz, k_rad, 8);
    #else
        SampleParticles(sim, sim.Lx, sim.Ly,         k_rad, 4);
    #endif
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    auto new_part_x = sim.particles.x;
    #ifdef THREEDIM
        TV tmp_part(sim.Lx/2.0, 0.0, 0.0);
    #else
        TV tmp_part(sim.Lx/2.0, 0.0);
    #endif
    new_part_x.push_back(tmp_part);
    sim.Np += 1;
    sim.particles = Particles(sim.Np);
    sim.particles.x = new_part_x;
    sim.delete_last_particle = false;

    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;

    sim.elastic_model = StvkWithHencky;
    sim.plastic_model = PerzynaDP;

    sim.use_von_mises_q = false;
    sim.use_pradhana = true;

    sim.perzyna_visc = 0;
    sim.perzyna_exp = 1;

    sim.dp_cohesion = 0;
    sim.dp_slope = std::tan(30.0 * M_PI / 180.0);

    std::string name;
    #ifdef THREEDIM
    name = "Back";    ObjectPlate sideback   = ObjectPlate(0,       1e20, -1e20, back,   SEPARATE, 0, name,   0,0,0,  1,0);  sim.plates.push_back(sideback);
    name = "Front";   ObjectPlate sidefront  = ObjectPlate(sim.Lz,  1e20, -1e20, front,  SEPARATE, 0, name,   0,0,0,  1,0);  sim.plates.push_back(sidefront);
    name = "Left";    ObjectPlate sideleft   = ObjectPlate(0,       1e20, -1e20, left,   SEPARATE, 0, name,   0,0,0,  1,0);  sim.plates.push_back(sideleft);
    name = "Ground";  ObjectPlate ground     = ObjectPlate(0,       1e20, -1e20, bottom, STICKY,   0, name,   0,0,0,  1,0);  sim.plates.push_back(ground);
    #else
    name = "Left";    ObjectPlate sideleft   = ObjectPlate(0,  1e20, -1e20, left,   SEPARATE, 0, name,   0,0,  1,0);  sim.plates.push_back(sideleft);
    name = "Ground";  ObjectPlate ground     = ObjectPlate(0,  1e20, -1e20, bottom, STICKY,   0, name,   0,0,  1,0);  sim.plates.push_back(ground);
    #endif

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x = (*max_x_it)(0);
    T diff = std::abs(max_x - 0.56);
    // 2D: 0.550019
    // 3D: 0.570681
    ASSERT_NEAR(diff, 0.0, 0.011);
}
