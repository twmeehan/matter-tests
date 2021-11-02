#include "simulation.hpp"

// TODO:
//  * CUDA Parallilize
//////////////////////////////////////////////////////////////////

// REMEMBER TO SPECIFY DIMENSION AND SPLINE DEGREE IN TOOLS!!!!!

int main(){

      Simulation sim;
      sim.sim_name = "dp_pradhana_xy_hard_novisc"; // "perzyna_dp_wieckowski_visc1";
      sim.directory = "/media/blatny/harddrive4/larsie/";
      sim.end_frame = 400;
      sim.fps = 100;

      T theta = 30 * M_PI / 180;
      sim.gravity = TV::Zero();
      sim.gravity[0] = +9.81 * std::sin(theta);
      sim.gravity[1] = -9.81 * std::cos(theta);

      sim.cfl = 0.6;
      sim.dt_max_coeff = 0.4;
      sim.flip_ratio = 0.99;
      sim.n_threads = 8;

      sim.initialize(/* E */ 3e6, /* nu */ 0.3, /* rho */ 2500);

      std::string sample = "samples/sample_Lx15.0_Ly3.0_r0.05/";
      std::ifstream file1(sample + "num.txt"); file1 >> sim.Np; file1.close();
      std::ifstream file2(sample + "Lx.txt");  file2 >> sim.Lx; file2.close();
      std::ifstream file3(sample + "Ly.txt");  file3 >> sim.Ly; file3.close();
      debug("Lx = ", sim.Lx);
      debug("Ly = ", sim.Ly);

      unsigned int Npx = std::sqrt(sim.Lx/sim.Ly * sim.Np);
      T dx_p = (sim.Lx / Npx);
      sim.dx = 2 * dx_p;

      sim.particles = Particles(sim.Np);
      unsigned int Np_check = load_array(sim.particles.x, sample + "pos.txt");
      if (Np_check != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }
      // for(int p = 0; p < sim.Np; p++){
      //     sim.particles.x[p][0] -= 1.0;
      // }

      sim.particle_volume = dx_p * dx_p;
      sim.particle_mass = sim.rho * sim.particle_volume;

      T vel_top   = 0.0;
      T vel_bot   = 0.0;
      T vel_left  = 0.0;
      T vel_right = 0.0;
     // T vel_back = 0.0;
     // T vel_front = 0.0;

      sim.vmin_factor = 25;
      sim.load_factor = 75;

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
      name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0 - offset,      0,                vel_left,  0,         1,               0,               left,   SEPARATE, name);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0,                vel_right, 0,         1,               0,               right,  SLIP, name);  sim.objects.push_back(side_right);

      sim.friction = 0.0; // in 3D currently only support zero friction, in 2D it's fine

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = NoPlasticity;
      // sim.plastic_model = VonMises;
      // sim.plastic_model = DruckerPrager;
      // sim.plastic_model = PerzynaVM;
      // sim.plastic_model = PerzynaDP;
      sim.plastic_model = PerzynaMuIDP;
      // sim.plastic_model = Curved;

      sim.dp_slope = 0.5;
      sim.dp_cohesion = 0;

      sim.yield_stress_orig = 1e3;
      sim.yield_stress_min = 1e3;

      sim.perzyna_exp = 1;
      sim.perzyna_visc = 0.0; // 50 s^-1

      sim.beta = 0.0;
      sim.M = 0.5;
      sim.p0 = 3461;

      sim.xi = 1;
      sim.xi_nonloc = 0;

      sim.nonlocal_l = 0;

      // For MCC only:
      // T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
      // sim.particles.eps_pl_vol_mcc.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_mcc.begin(), sim.particles.eps_pl_vol_mcc.end(), eps_pl_vol_init );

      sim.simulate();
      // sim.validateRMA();


	return 0;
}
