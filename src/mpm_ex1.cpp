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
      sim.sim_name = "reg_test";
      sim.end_frame = 100;
      sim.frame_dt = 1.0 / 300.0;
      sim.gravity = TV::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.4;
      sim.flip_ratio = 0.0;
      sim.n_threads = 24;

      sim.Lx = 1.0;
      sim.Ly = 2.0;
      sim.Lz = 1.0;

      int Npx = 51;
      int Npy = 101;
      int Npz = 51;

      T dxp = sim.Lx / (Npx-1.0);

      sim.dx = 2.0 * dxp;

      const unsigned int ppc = 8;  // 2D: ppc = 4
      sim.Np = Npx * Npy * Npz;

      // sim.Np = 1762;
      // sim.Np = 154360;
      // sim.Np = 507022;

      // Convention below: The conditional operator returns TRUE if DILATION
      T vel_top = -0.2;
      T vel_bot = 0.2;
      std::string name;
      name = "Compressor"; InfinitePlate compressor = InfinitePlate(0, sim.Ly, 0, 0, vel_top, 0,    top,  STICKY, name);  sim.objects.push_back(compressor);
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0, 0,      0, 0, vel_bot, 0, bottom,  STICKY, name);  sim.objects.push_back(ground);
      sim.friction = 0.0; // currently only support zero friction

      sim.initialize(/* E */ 5e7, /* nu */ 0.0, /* rho */ 1000);
      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass   = ", sim.particle_mass);
      debug("Np              = ", sim.Np);

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      sim.plastic_model = VonMises;
      sim.xi = 0; // for both VM and DP

      // Von Mises:
      sim.yield_stress_orig = std::sqrt(2.0/3.0) * /* q_max */ 2e6;
      sim.yield_stress_min  = std::sqrt(2.0/3.0) * /* q_max */ 2e6;

      // DPSimpleSoft
      sim.friction_angle = 13.342363799593;
      sim.cohesion = 0.017e6 / (sim.K * sim.dim); // p_min = - K * dim * cohesion

      // Regularization by Laplacian
      sim.reg_length = 0; //sim.dx;

      // Random samples from file
      //unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/phd-stuff/gold/output/microstructures/benchmarks/v4_N10000_phi03_mesh40_mc6_Lrve1_b50_mu1_typev_seed42/xyz_ext8_rand6.txt");
      // unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/larsiempm/build/microstructures/m65_mc9_phi026_seed12/xyz.txt");
      // debug("Np (load_array)  = ", Np);
      // if (Np != sim.Np){
      //     debug("Particle number mismatch!!!");
      //     return 0;
      // }

      // Initial state
      sim.amplitude = 0.0;

      // 2D:
      // std::vector<T> disp_i(ppc); disp_i = {0.25, 0.75, 0.25, 0.75};
      // std::vector<T> disp_j(ppc); disp_j = {0.25, 0.75, 0.75, 0.25};

      // 3D
      std::vector<T> disp_i(ppc); disp_i = {0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75};
      std::vector<T> disp_j(ppc); disp_j = {0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25};
      std::vector<T> disp_k(ppc); disp_k = {0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75};

      int p = -1;
      for(int i = 0; i < Npx; i++){
          for(int j = 0; j < Npy; j++){
              for(int k = 0; k < Npz; k++){
                  p++;

                  // T px = (i+disp_i[d])*sim.dx;
                  // T py = (j+disp_j[d])*sim.dx;
                  // T pz = (k+disp_k[d])*sim.dx; // 2D: pz = 0

                  T px = (i)*dxp;
                  T py = (j)*dxp;
                  T pz = (k)*dxp; // 2D: pz = 0

                  // CASE 3:
                  // T pvx = sim.amplitude*std::sin( M_PI*(px-0.5) );
                  // T pvy = sim.amplitude*std::sin( M_PI*(py-0.5) );
                  // CASE 2:
                  // T pvx = sim.amplitude * sim.frame_dt * sim.end_frame * px;
                  // T pvy = sim.amplitude * sim.frame_dt * sim.end_frame * py ;
                  // CASE 1:
                  // T pvx = sim.amplitude * sim.frame_dt * sim.end_frame;
                  // T pvy = sim.amplitude * sim.frame_dt * sim.end_frame;
                  // CASE 0:
                  T pvx = sim.amplitude;
                  T pvy = sim.amplitude;
                  T pvz = sim.amplitude;
                  sim.particles.x[p](0) = px;
                  sim.particles.x[p](1) = py;
                  sim.particles.x[p](2) = pz;
                  sim.particles.v[p](0) = pvx;
                  sim.particles.v[p](1) = pvy;
                  sim.particles.v[p](2) = pvz;
              } // end for k
          } // end for i
      } // end for j
      debug("Added particles = ", p);
      if ((p+1) != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }

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
