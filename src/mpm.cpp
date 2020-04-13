#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * 3D
//  * Complete doxygen
//  * Alembic output
//  * Boundary condition for Laplace? Basically using Dirichlet now
//  * OMP-Parallilize explicit_euler_update
//  * CUDA Parallilize
//  * Position and def. grad update should be based on velocity before application of friction (but after collision). Problem for SLIP, not STICKY
//////////////////////////////////////////////////////////////////

int main(){

      Simulation sim;
      sim.sim_name = "threedim-micro-plastic-2";
      sim.end_frame = 80;
      sim.frame_dt = 1.0 / 800.0;
      sim.dx = 1.0 / 40;
      sim.L = 1;
      sim.gravity = TV::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;
      sim.flip_ratio = 0.9;
      sim.n_threads = 24;

      // const unsigned int Nloop = std::round(1.0/sim.dx);
      // const unsigned int ppc  = 8;          // 2D: ppc = 4
      // sim.Np = Nloop * Nloop * Nloop * ppc; // 2D: Np = Nloop * Nloop * ppc;
      // sim.Np = 1762;
      // sim.Np = 154360;
      sim.Np = 19295;

      std::string name;
      name = "Ground";     InfinitePlate ground      = InfinitePlate(0, 0, 0, 0, +2,  0, bottom, SLIP, name);  sim.objects.push_back(ground);
      name = "Compressor"; InfinitePlate compressor  = InfinitePlate(0, 1, 0, 0, -2, 0,  top,    SLIP, name);  sim.objects.push_back(compressor);
      name = "Left";       InfinitePlate left_plate  = InfinitePlate(0, 0, 0, +2, 0,  0, left,   SLIP, name);  sim.objects.push_back(left_plate);
      name = "Right";      InfinitePlate right_plate = InfinitePlate(1, 0, 0, -2, 0,  0, right,  SLIP, name);  sim.objects.push_back(right_plate);
      name = "Front";      InfinitePlate front_plate = InfinitePlate(0, 0, 0, 0, 0,  +2, front,  SLIP, name);  sim.objects.push_back(front_plate);
      name = "Back";       InfinitePlate back_plate  = InfinitePlate(0, 0, 1, 0, 0,  -2, back,   SLIP, name);  sim.objects.push_back(back_plate);

      sim.friction = 0.0;

      sim.initialize(/* E */ 1e8, /* nu */ 0.3, /* rho */ 1000*10);
      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass   = ", sim.particle_mass);
      debug("Np              = ", sim.Np);

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      sim.plastic_model = DPSimpleSoft;
      sim.xi = 1e15; // for both VM and DP

      // Von Mises:
      sim.yield_stress_orig = std::sqrt(2.0/3.0) * /* q_max */ 4000.0;
      sim.yield_stress_min  = std::sqrt(2.0/3.0) * /* q_max */ 4000.0;

      // DPSimpleSoft
      sim.friction_angle = 13.34;
      sim.cohesion = 0.017e6 / (sim.K * sim.dim); // p_min = - K * dim * cohesion

      // Regularization by Laplacian
      sim.reg_length = sim.dx;
      sim.reg_const = 0; // 2*sim.mu;

      // Random samples from file
      //unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/phd-stuff/gold/output/microstructures/benchmarks/v4_N10000_phi03_mesh40_mc6_Lrve1_b50_mu1_typev_seed42/xyz_ext8_rand6.txt");
      unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/phd-stuff/gold/output/microstructures/benchmarks/v4_N10000_phi03_mesh40_mc6_Lrve1_b50_mu1_typev_seed42/xyz.txt");
      debug("Np load_array)  = ", Np);
      if (Np != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }

      // // Initial state
      // sim.amplitude = 0.0;

      // 2D:
      // std::vector<T> disp_i(ppc); disp_i = {0.25, 0.75, 0.25, 0.75};
      // std::vector<T> disp_j(ppc); disp_j = {0.25, 0.75, 0.75, 0.25};

      // 3D
      // std::vector<T> disp_i(ppc); disp_i = {0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75};
      // std::vector<T> disp_j(ppc); disp_j = {0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25};
      // std::vector<T> disp_k(ppc); disp_k = {0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75};
      //
      // int p = -1;
      // for(int i = 0; i < Nloop; i++){
      //     for(int j = 0; j < Nloop; j++){
      //         for(int k = 0; k < Nloop; k++){ // 2D: Nloop must be 1
      //             for(int d = 0; d < ppc; d++){
      //                 p++;
      //                 T px = (i+disp_i[d])*sim.dx;
      //                 T py = (j+disp_j[d])*sim.dx;
      //                 T pz = (k+disp_k[d])*sim.dx; // 2D: pz = 0
      //                 // CASE 3:
      //                 // T pvx = sim.amplitude*std::sin( M_PI*(px-0.5) );
      //                 // T pvy = sim.amplitude*std::sin( M_PI*(py-0.5) );
      //                 // CASE 2:
      //                 // T pvx = sim.amplitude * sim.frame_dt * sim.end_frame * px;
      //                 // T pvy = sim.amplitude * sim.frame_dt * sim.end_frame * py ;
      //                 // CASE 1:
      //                 // T pvx = sim.amplitude * sim.frame_dt * sim.end_frame;
      //                 // T pvy = sim.amplitude * sim.frame_dt * sim.end_frame;
      //                 // CASE 0:
      //                 T pvx = sim.amplitude;
      //                 T pvy = sim.amplitude;
      //                 T pvz = sim.amplitude;
      //                 sim.particles.x[p](0) = px;
      //                 sim.particles.x[p](1) = py;
      //                 sim.particles.x[p](2) = pz;
      //                 sim.particles.v[p](0) = pvx;
      //                 sim.particles.v[p](1) = pvy;
      //                 sim.particles.v[p](2) = pvz;
      //             } // end for d
      //         } // end for k
      //     } // end for i
      // } // end for j
      // debug("Added particles = ", p);
      // if ((p+1) != sim.Np){
      //     debug("Particle number mismatch!!!");
      //     return 0;
      // }

    sim.simulate();


    ///////////// ALEMBIC TESTING: ////////////////
    /*
    csv2abc("test.abc");
    return 0;
    */
    ///////////// DEBUG 1: ////////////////
    /*
        sim.P2G();
        sim.saveSim();
        sim.saveGridVelocities();
        sim.current_time_step++;

        for (int i=0; i<50; i++){
            sim.G2P();
            sim.P2G();
            sim.saveSim();
            sim.saveGridVelocities();
            sim.current_time_step++;
        }
    */
    ///////////// DEBUG 2: ////////////////
    /*
        for (int i=0; i<50; i++){
            sim.updateDt();
            sim.P2G();
            sim.saveGridVelocities("before_");
            sim.explicitEulerUpdate();
            sim.saveGridVelocities("after_");
            sim.G2P();
            sim.saveSim("before_");
            sim.deformationUpdate();
            sim.positionUpdate();
            sim.saveSim("after_");
            sim.current_time_step++;
        }
    */


	return 0;
}
