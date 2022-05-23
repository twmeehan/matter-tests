#include "simulation.hpp"

// TODO:
//  * CUDA Parallilize

// REMEMBER TO SPECIFY DIMENSION AND SPLINE DEGREE IN TOOLS!!!!!

int main(){

      Simulation sim;
      sim.sim_name = "sintering_test_compress";
      // sim.sim_name = "dpmui_Lx1.6";
      sim.directory = "/home/multiples/Documents/larsie/";
      sim.end_frame = 500; // 600 // 400
      sim.fps = 20; // 100 // 4

      T theta = 0 * M_PI / 180;
      sim.gravity = TV::Zero();
      // sim.gravity[0] = +9.81 * std::sin(theta);
      // sim.gravity[1] = -9.81 * std::cos(theta);
      sim.gravity_time = 0.0;

      sim.cfl = 0.6;
      sim.dt_max_coeff = 0.45;
      sim.flip_ratio = 0.99; //0.99;
      sim.n_threads = 8;

      // sim.initialize(/* E */ 1e9, /* nu */ 0.3, /* rho */ 1630);
      // sim.initialize(/* E */ 1e5, /* nu */ 0.3, /* rho */ 1500);
      sim.initialize(/* E */ 1e5, /* nu */ 0.3, /* rho */ 370);

      sim.Lx = 0.25;
      sim.Ly = 0.2;
      // sim.Lz = 0.35;

      SampleInBox(sim.Lx, sim.Ly,         0.0012, sim);
      // SampleInBox(sim.Lx, sim.Ly, sim.Lz, 0.01, sim);

      T offset = -0.01 * sim.dx/2.0; // When the grid is aligned with the boundary, it is important that the object overlap a bit into the particle domain
      std::string name;
      ///////// 3D /////////
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0-offset,     1e20, -1e20, bottom, STICKY,   0.51, name,   0,0,0,   1,0);  sim.objects.push_back(ground);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0-offset,     1e20, -1e20, left,   SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx+offset, 0.2,   0,   right,  SEPARATE, 0.18, name,   0,3,0,   1,0);  sim.objects.push_back(side_right);
      // name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0-offset,     1e20, -1e20, back,   SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_back);
      // name = "SideFront";  InfinitePlate side_front = InfinitePlate(sim.Lz+offset,1e20, -1e20, front,  SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_front);
      ///////// 2D /////////
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0-offset,      1e20, -1e20, bottom, SLIP, 0.5, name,   0, 0,     1,0);  sim.objects.push_back(ground);
      name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0-offset,      1e20, -1e20, left,   SLIP, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_left);
      name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx+offset, 1e20, -1e20, right,  SLIP, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_right);
      name = "Compressor"; InfinitePlate compressor = InfinitePlate(sim.Ly+offset, 1e20, -1e20, top,    SEPARATE, 0.0, name,   0, -0.005,  1,0);  sim.objects.push_back(compressor);

      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = NoPlasticity;
      // sim.plastic_model = VonMises;
      // sim.plastic_model = DruckerPrager;
      // sim.plastic_model = PerzynaVM;
      // sim.plastic_model = PerzynaDP;
      // sim.plastic_model = PerzynaMuIDP;
      sim.plastic_model = ModifiedCamClay;
      // sim.plastic_model = ModifiedCamClayHard;
      // sim.plastic_model = PerzynaMCC;
      // sim.plastic_model = PerzynaMuIMCC;
      // sim.plastic_model = PerzynaSinterMCC;

      sim.dp_slope = 0.5;
      sim.dp_cohesion = 0;

      sim.vm_ptensile = -5e10;
      sim.yield_stress_orig = 5e3;
      sim.yield_stress_min = sim.yield_stress_orig;

      sim.perzyna_exp = 1;
      sim.perzyna_visc = 0.001;

      sim.beta = 0.2;
      sim.M = sim.dp_slope; //0.6;
      sim.p0 = 1e2;

      sim.xi = 0.05;  // In case of sintering, xi = 1/(lambda*phi_0)
      sim.xi_nonloc = 0;

      sim.nonlocal_l = 0;

      // For mu(I) rheology only
      sim.rho_s           = 2450; // 2672.13;
      sim.grain_diameter  = 7e-4;
      sim.in_numb_ref     = 1e-3; // 3.892e-2;
      sim.mu_1            = sim.dp_slope; //sim.M
      sim.mu_2            = 0.7;

      // For sintering only
      sim.sinter_tc = 20;   // 1     // 20
      sim.sinter_ec = 0.03; // 1e-3  // 0.03
      sim.sinter_Sinf = 20;
      sim.particles.sinter_S.resize(sim.Np); std::fill( sim.particles.sinter_S.begin(), sim.particles.sinter_S.end(), sim.sinter_Sinf );

      // For MCC only:
      T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
      sim.particles.eps_pl_vol_mcc.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_mcc.begin(), sim.particles.eps_pl_vol_mcc.end(), eps_pl_vol_init );

      sim.simulate();
      // sim.validateRMA();


	return 0;
}
