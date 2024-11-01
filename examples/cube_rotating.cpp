#include "simulation.hpp"
#include "tools.hpp"
#include "sampling_particles.hpp"

int main(){

    Simulation sim;

    sim.directory = "output/";
    sim.sim_name = "cube_rotating";
    sim.end_frame = 200;
    sim.fps = 10;

    sim.save_grid = true;

    sim.gravity = TV::Zero();

    sim.cfl = 0.5;
    sim.flip_ratio = -1;
    sim.n_threads = 8;

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1550);

    T h_gate, l_gate;
    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.05;
        SampleParticles(sim, sim.Lx, sim.Ly, sim.Lz, 0.01);
    #else
        SampleParticles(sim, sim.Lx, sim.Ly, 0.01);
    #endif

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

    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;

    // Elastoplasticity
    sim.elastic_model = StvkWithHencky;
    sim.plastic_model = NoPlasticity;

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
