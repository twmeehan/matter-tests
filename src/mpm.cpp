#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * Complete doxygen
//  * CUDA Parallilize
//  * Position and def. grad update should be based on velocity before application of friction (but after collision). Problem for SLIP, not STICKY
//////////////////////////////////////////////////////////////////

// REMEMBER TO SPECIFY DIMENSION IN TOOLS!!!!!

int main(){

      Simulation sim;
      sim.sim_name = "test_larsie_2"; // "perzyna_dp_wieckowski_visc1";
      sim.directory = "/media/blatny/harddrive4/larsie/"; // "dumps/";
      sim.end_frame = 500;
      sim.fps = 100;
      sim.gravity = TV::Zero(); sim.gravity[1] = -9.81;
      sim.cfl = 0.6;
      sim.dt_max_coeff = 0.4;
      sim.flip_ratio = 0.99;
      sim.n_threads = 16;

      sim.Lx = 2.0;
      sim.Ly = 2.0;
      // sim.Lz = sim.Lx;

      ///////////// PARTICLES FROM FILE //////////////////
      std::string sample = "samples/sample_w2_h2_r0.011_";
      std::ifstream file(sample + "num.txt");
      file >> sim.Np;
      file.close();
      unsigned int Npx = std::sqrt(sim.Lx/sim.Ly * sim.Np);
      sim.dx = 2 * (sim.Lx / Npx);

      ///////////// PARTICLES ON GRID //////////////////
      // int Npx = 30+1;
      // int Npy = 60+1;
      //
      // T dxp = sim.Lx / (Npx-1.0);
      // sim.dx = 2.0 * dxp;
      // sim.Np = Npx * Npy;// * Npz;
      ///////////////////////////////////////////////////

      T vel_top   = 0.0;
      T vel_bot   = 0.0;
      T vel_left  = 0.0;
      T vel_right = 0.0;
     // T vel_back = 0.0;
     // T vel_front = 0.0;

      sim.vmin_factor = 25;
      sim.load_factor = 75; //75;

      T offset = -0.01 * sim.dx/2.0; // When the grid is aligned with the boundary, it is important that the object overlap a bit into the particle domain

      std::string name;
      ///////// 3D /////////
      // name = "Compressor"; InfinitePlate compressor = InfinitePlate(0, sim.Ly + offset, 0,    0, vel_top, 0, sim.vmin_factor, sim.load_factor,    top, SLIP, name);  sim.objects.push_back(compressor);
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0, 0      - offset, 0,    0, vel_bot, 0,               1,               0, bottom, SLIP, name);  sim.objects.push_back(ground);
      //
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0      - offset, 0, 0,    vel_left,  0, 0,             1,               0,   left, SLIP, name);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0, 0,    vel_right, 0, 0,             1,               0,  right, SLIP, name);  sim.objects.push_back(side_right);
      //
      // name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0, 0, 0      - offset,    0, 0, vel_back,              1,               0,  back,  SLIP, name);  sim.objects.push_back(side_back);
      // name = "SideFront";  InfinitePlate side_front = InfinitePlate(0, 0, sim.Lz + offset,    0, 0, vel_front,             1,               0,  front, SLIP, name);  sim.objects.push_back(side_front);

      ///////// 2D /////////
      // name = "Compressor"; InfinitePlate compressor = InfinitePlate(0,               sim.Ly + offset,  0,         vel_top,   sim.vmin_factor, sim.load_factor, top,    SLIP, name);  sim.objects.push_back(compressor);
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0,               0 - offset,       0,         vel_bot,   1,               0,               bottom, STICKY, name);  sim.objects.push_back(ground);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0 - offset,      0,                vel_left,  0,         1,               0,               left,   SLIP, name);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0,                vel_right, 0,         1,               0,               right,  SLIP, name);  sim.objects.push_back(side_right);

      sim.friction = 0.0; // currently only support zero friction

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = NoPlasticity;
      // sim.plastic_model = VonMises;
      // sim.plastic_model = DruckerPrager;
      // sim.plastic_model = PerzynaVM;
      sim.plastic_model = PerzynaNA;
      // sim.plastic_model = Curved;

      sim.dp_slope = 1.02857; // 30 deg
      sim.dp_cohesion = 0;

      sim.yield_stress_orig = 1e2;
      sim.yield_stress_min = 1e2;

      sim.perzyna_exp = 1;
      sim.perzyna_visc = 1.0 / 50; // 50 s^-1

      sim.beta = 0.2;
      sim.M = 0.1;
      sim.p0 = 1e3;

      sim.xi = 0;
      sim.xi_nonloc = 0;

      sim.nonlocal_l = 0;

      sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1500);
      debug("Wave speed       = ", sim.wave_speed);
      debug("dt_max           = ", sim.dt_max);
      debug("particle_volume  = ", sim.particle_volume);
      debug("particle_mass    = ", sim.particle_mass);
      debug("Np               = ", sim.Np);

      // T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
      // sim.particles.eps_pl_vol_3.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_3.begin(), sim.particles.eps_pl_vol_3.end(), eps_pl_vol_init );

      ////////////////////////////////////////////////////
      ///////////// PARTICLES FROM FILE //////////////////
      ////////////////////////////////////////////////////
      unsigned int Np_check = load_array(sim.particles.x, sample + "pos.txt");
      if (Np_check != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }

      for(int p = 0; p < sim.Np; p++){
          sim.particles.x[p][0] -= 1.0;
      }

      //////////////////////////////////////////////////
      ///////////// PARTICLES ON GRID //////////////////
      //////////////////////////////////////////////////
      // int p = -1;
      // for(int i = 0; i < Npx; i++){
      //     for(int j = 0; j < Npy; j++){
      //       //  for(int k = 0; k < Npz; k++){
      //             p++;
      //
      //             // T px = (i+disp_i[d])*sim.dx;
      //             // T py = (j+disp_j[d])*sim.dx;
      //             T px = (i)*dxp;
      //             T py = (j)*dxp;
      //            // T pz = (k)*dxp;
      //
      //             sim.particles.x[p](0) = px;
      //             sim.particles.x[p](1) = py;
      //           // sim.particles.x[p](2) = pz;
      //
      //             sim.particles.v[p](0) = 0;
      //             sim.particles.v[p](1) = 0;
      //            // sim.particles.v[p](2) = 0;
      //
      //             // sim.particles.yield_stress_orig[p] = ys;
      //
      //       //  } // end for k
      //     } // end for j
      // } // end for i
      // debug("Added particles = ", p);
      // if ((p+1) != sim.Np){
      //     debug("Particle number mismatch!!!");
      //     return 0;
      // }
      /////////////////////////////////////////////////////

      sim.simulate();
      // sim.validateRMA();


	return 0;
}
