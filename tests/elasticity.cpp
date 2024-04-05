#include <gtest/gtest.h>

#include "../src/simulation.hpp"
#include "../src/tools.hpp"
#include "../src/sampling/sampling_particles.hpp"

// Demonstrate some basic assertions.
TEST(MyTest, BasicAssertions) {

    Simulation sim;
    sim.save_sim = false;
    sim.end_frame = 25;
    sim.fps = 100;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 0;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 2000); // feeder, incl slope and gran collapse


    sim.Lx = 0.1;
    sim.Ly = 0.1;
    T k_rad = 0.001;
    SampleParticles(sim.Lx, sim.Ly,     k_rad, 6, 0, sim);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) -= 0.5*sim.Ly;
    }

    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;
    T vel = 0.1;

    std::string name;
    name = "Ground";  ObjectPlate ground = ObjectPlate(-0.5*sim.Lx - 0.5*sim.dx,  1e20, -1e20, bottom, STICKY, 0, name,   0,  vel,     1,0);  sim.plates.push_back(ground);
    name = "Compre";  ObjectPlate compre = ObjectPlate( 0.5*sim.Ly + 0.5*sim.dx,  1e20, -1e20, top,    STICKY, 0, name,   0, -vel,     1,0);  sim.plates.push_back(compre);

    sim.elastic_model = StvkWithHencky;
    sim.plastic_model = NoPlasticity;

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    auto min_x_it = std::min_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0); } );
    T max_x = (*max_x_it)(0);
    T min_x = (*min_x_it)(0);

    debug(max_x);
    debug(min_x);

  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}
