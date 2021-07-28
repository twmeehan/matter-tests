#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * Complete doxygen
//  * CUDA Parallilize
//  * Position and def. grad update should be based on velocity before application of friction (but after collision). Problem for SLIP, not STICKY
//////////////////////////////////////////////////////////////////

int main(){

      Simulation sim;
      sim.sim_name = "3d_ql_anal_soft_rho300_xi0.3";
      // sim.sim_name = "3d_elastic";
      // sim.sim_name = "test_rma_quad_anal";
      sim.directory = "/media/blatny/harddrive4/larsie/"; // "dumps/";
      sim.end_frame = 800;
      sim.fps = 120;
      sim.gravity = TV::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;
      sim.dt_max_coeff = 0.4;
      sim.flip_ratio = 0.99;
      sim.n_threads = 24;

      sim.Lx = 1.0;
      sim.Ly = 2.0;
      sim.Lz = sim.Lx;

      int Npx = 30+1;
      int Npy = 60+1;
      int Npz = Npx;

      T dxp = sim.Lx / (Npx-1.0);

      sim.dx = 2.0 * dxp;

      sim.Np = Npx * Npy * Npz;

      T vel_top = -0.2;
      T vel_bot = 0.0;
      T vel_left = 0.0;
      T vel_right = 0.0;
      T vel_back = 0.0;
      T vel_front = 0.0;

      T vmin_factor = 25;
      T load_factor = 75;

      T offset = 0.1 * dxp;

      std::string name;
      name = "Compressor"; InfinitePlate compressor = InfinitePlate(0, sim.Ly + offset, 0,    0, vel_top, 0,     vmin_factor, load_factor,    top, SLIP, name);  sim.objects.push_back(compressor);
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0, 0      - offset, 0,    0, vel_bot, 0,               1,           0, bottom, SLIP, name);  sim.objects.push_back(ground);

      name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0      - offset, 0, 0,    vel_left,  0, 0,             1,           0,   left, SLIP, name);  sim.objects.push_back(side_left);
      name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0, 0,    vel_right, 0, 0,             1,           0,  right, SLIP, name);  sim.objects.push_back(side_right);

      name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0, 0, 0      - offset,    0, 0, vel_back,              1,           0,  back,  SLIP, name);  sim.objects.push_back(side_back);
      name = "SideFront";  InfinitePlate side_front = InfinitePlate(0, 0, sim.Lz + offset,    0, 0, vel_front,             1,           0,  front, SLIP, name);  sim.objects.push_back(side_front);

      sim.friction = 0.0; // currently only support zero friction

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = VonMises;
      sim.plastic_model = QuadraticLars;
      // sim.plastic_model = NoPlasticity;
      sim.beta = 0.43-0.40*0.3;
      sim.M = 1.35;
      sim.p0 = 50e3;

      sim.xi = 0.3; // 1e20;

      sim.xi_nonloc = 0;
      sim.nonlocal_l = 0;

      T ys = 50e3;



      sim.initialize(/* E */ 3e8, /* nu */ 0.3, /* rho */ 300);
      debug("Wave speed       = ", sim.wave_speed);
      debug("dt_max           = ", sim.dt_max);
      debug("particle_volume  = ", sim.particle_volume);
      debug("particle_mass    = ", sim.particle_mass);
      debug("Np               = ", sim.Np);

      // Random samples from file
      //unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/phd-stuff/gold/output/microstructures/benchmarks/v4_N10000_phi03_mesh40_mc6_Lrve1_b50_mu1_typev_seed42/xyz_ext8_rand6.txt");
      // unsigned int Np = load_array(sim.particles.x, "/home/blatny/repos/larsiempm/build/microstructures/m65_mc9_phi026_seed12/xyz.txt");
      // debug("Np (load_array)  = ", Np);
      // if (Np != sim.Np){
      //     debug("Particle number mismatch!!!");
      //     return 0;
      // }

      int p = -1;
      for(int i = 0; i < Npx; i++){
          for(int j = 0; j < Npy; j++){
              for(int k = 0; k < Npz; k++){
                  p++;

                  // T px = (i+disp_i[d])*sim.dx;
                  // T py = (j+disp_j[d])*sim.dx;
                  T px = (i)*dxp;
                  T py = (j)*dxp;
                  T pz = (k)*dxp;

                  sim.particles.x[p](0) = px;
                  sim.particles.x[p](1) = py;
                  sim.particles.x[p](2) = pz;

                  sim.particles.v[p](0) = 0;
                  sim.particles.v[p](1) = 0;
                  sim.particles.v[p](2) = 0;

                  sim.particles.yield_stress_orig[p] = ys;

              } // end for k
          } // end for j
      } // end for i
      debug("Added particles = ", p);
      if ((p+1) != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }

    sim.simulate();
    // sim.validateRMA();

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
