#include "simulation.hpp"
#include "sampling_particles.hpp"

// TODO:
//  * CUDA Parallilize

// REMEMBER TO SPECIFY DIMENSION AND SPLINE DEGREE IN TOOLS!!!!!

int main(){

      Simulation sim;

      sim.directory = "/media/blatny/harddrive4/larsie/";
      sim.end_frame = 200;
      // sim.fps = 238.2166843; // gran collapse, so t_star at frame 100
      // sim.fps = 50; // inc slope
      sim.fps = 7; // 25

      // sim.sim_name = "pbc_a28_dx4_jop_E1e9_d7mm";
      sim.sim_name = "feeder_kam_a37_d2-5mm_dx35_Lx18";
      // sim.sim_name = "rounded_dp_kam_a22_d500mm";
      // sim.sim_name = "inclslope0_mu25_kam_mccmui";
      // sim.sim_name = "grancoll_kam_mccmui_vgate045";

      T theta_deg = 37;
      T theta = theta_deg * M_PI / 180;
      sim.gravity = TV::Zero();
      sim.gravity[0] = +9.81 * std::sin(theta);
      sim.gravity[1] = -9.81 * std::cos(theta);
      sim.gravity_time = 0.0;

      sim.cfl = 0.5;
      sim.flip_ratio = -0.95;
      sim.n_threads = 2;

      // sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1450); // gran collapse
      sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1550); // incl slope and feeder

      // sim.pbc = true;
      // sim.Lx = 12*0.00125;
      // sim.Ly = 0.05;
      // SampleParticles(sim.Lx, sim.Ly, 0.00023, 20, 0, sim);
      // sim.dx = 0.00125;
      // sim.particle_volume = sim.Lx * sim.Ly / sim.Np;
      // sim.particle_mass = sim.rho * sim.particle_volume;

      // sim.pbc = true;
      // sim.Np = 20;
      // sim.Ly = 1;
      // sim.dx = sim.Ly / sim.Np; // one particle per cell
      // sim.Lx = sim.dx;
      // sim.particle_volume = sim.dx*sim.dx;
      // sim.particle_mass = sim.rho * sim.particle_volume;
      // sim.particles = Particles(sim.Np);
      // for(int p = 0; p < sim.Np; p++){
      //     TV vec = TV::Ones()*0.5*sim.dx;
      //     vec(1) = (0.5+p)*sim.dx;
      //     sim.particles.x[p] = vec;
      // }

      ///// Incl slope for h_stop (Rounded Edges)
      // sim.pbc = false;
      // sim.Lx = 0.9;
      // sim.Ly = 0.1;
      // SampleParticles(sim.Lx, sim.Ly,     0.001, 6, 4, sim);

      ///// Gran collapse:
      // sim.pbc = false;
      // sim.Lx = 0.08;
      // sim.Ly = 0.05;
      // SampleParticles(sim.Lx, sim.Ly,     0.0004, 6, 0, sim);

      ///// Incl slope:
      // sim.pbc = false;
      // sim.Lx = 0.20;
      // sim.Ly = 0.14;
      // SampleParticles(sim.Lx, sim.Ly,     0.001, 6, 0, sim);

      ///// Feeder with ramp
      sim.pbc = false;
      sim.Lx = 1.8;
      sim.Ly = 0.1;
      T h_gate = 0.05;
      T l_gate = 0.1;
      // SampleParticles(sim.Lx, sim.Ly,     0.002, 6, 0, sim);
      // SampleParticles(sim.Lx, sim.Ly,     0.0005, 6, 0, sim); // dx3
      SampleParticles(sim.Lx, sim.Ly,     0.0004, 6, 0, sim); // dx35
      // SampleParticles(sim.Lx, sim.Ly,     0.0003, 6, 0, sim); // dx4
      for(int p = 0; p < sim.Np; p++){
          sim.particles.x[p](0) -= sim.Lx;
          sim.particles.x[p](1) += l_gate + 0.5*sim.dx;
      }
      auto new_part_x = sim.particles.x;
      TV tmp_part(l_gate, 0.0);
      new_part_x.push_back(tmp_part);
      sim.Np += 1;
      sim.particles = Particles(sim.Np);
      sim.particles.x = new_part_x;


      sim.dt_max = 0.5 * sim.dx / sim.wave_speed;
      // sim.dt_max = 0.5 * sim.dx / (std::sqrt(1.0e7/sim.rho));

      T offset = 0; //-0.1 * sim.dx/2.0; // When the grid is aligned with the boundary, it is important that the object overlap a bit into the particle domain
      std::string name;
      ///////// 3D /////////
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0-offset,     1e20, -1e20, bottom, STICKY,   0.51, name,   0,0,0,   1,0);  sim.objects.push_back(ground);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0-offset,     1e20, -1e20, left,   SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx+offset, 0.2,   0,   right,  SEPARATE, 0.18, name,   0,3,0,   1,0);  sim.objects.push_back(side_right);
      // name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0-offset,     1e20, -1e20, back,   SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_back);
      // name = "SideFront";  InfinitePlate side_front = InfinitePlate(sim.Lz+offset,1e20, -1e20, front,  SEPARATE, 0.18, name,   0,0,0,   1,0);  sim.objects.push_back(side_front);
      ///////// 2D /////////
      //// Granular collapse
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0-offset,1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,        1,0);  sim.objects.push_back(ground);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,       1e20, 0,     left,   SEPARATE, 0.0, name,   0, 0.45,     1,0);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx,  1e20, 0,     right,  SEPARATE, 0.0, name,   0, 0.45,     1,0);  sim.objects.push_back(side_right);
      //// Inclined slope
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0-offset,1e20, -1e20, bottom, STICKY,   0.0, name,  0, 0,    1,0);  sim.objects.push_back(ground);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      1e20, -1e20,  left,   SEPARATE, 0.0, name,  0, 0,    1,0);  sim.objects.push_back(side_left);
      //// Feeder
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0,      1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,     1,0);  sim.objects.push_back(ground);
      name = "SideRight";  InfinitePlate side_right = InfinitePlate(l_gate, 1e20, h_gate, right, SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_right);
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      1e20, -1e20, left,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_left);
      // name = "Top";        InfinitePlate lid        = InfinitePlate(h_gate, 1e20, l_gate, top,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(lid);
      name = "Ramp";       InfinitePlate ramp       = InfinitePlate(l_gate,      0, -1e20, bottom, SEPARATE, 0.6, name,   0, 0,     1,0);  sim.objects.push_back(ramp);
      name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      l_gate, -1e20, left,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_left);


      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = NoPlasticity;
      // sim.plastic_model = VonMises;
      // sim.plastic_model = DruckerPrager;
      // sim.plastic_model = PerzynaVM;
      // sim.plastic_model = PerzynaDP;
      // sim.plastic_model = PerzynaMuIDP;
      // sim.plastic_model = MCC;
      // sim.plastic_model = MCCHard;
      // sim.plastic_model = MCCHardExp;
      // sim.plastic_model = PerzynaMCC;
      // sim.plastic_model = PerzynaMCCHard;
      sim.plastic_model = PerzynaMuIMCC;
      // sim.plastic_model = PerzynaSinterMCC;

      sim.use_jop_definitions = true;

      // sim.dp_slope = 3*std::sqrt(3.0) / 3.0 * std::tan(20.9 * M_PI / 180);
      sim.dp_slope = std::tan(20.9 * M_PI / 180);
      // sim.dp_slope = std::tan(25 * M_PI / 180); // = angle of repose in incl slope test
      sim.dp_cohesion = 0;

      sim.vm_ptensile = -5e10;
      sim.yield_stress_orig = 5e3;
      sim.yield_stress_min = sim.yield_stress_orig;

      sim.perzyna_exp = 1;
      sim.perzyna_visc = 1;

      sim.beta = 0.0;
      sim.M = sim.dp_slope;
      sim.p0 = 100;

      sim.xi = 50; // In case of sintering, xi = 1/(lambda*phi_0)

      sim.xi_nonloc = 0;
      sim.nonlocal_l = 0;

      // For mu(I) rheology only
      sim.rho_s           = 2500;
      sim.grain_diameter  = 2.5e-3; // incl slope and feeder
      // sim.grain_diameter  = 2e-3; // gran collapse
      sim.in_numb_ref     = 0.279;
      sim.mu_1            = sim.dp_slope; //sim.M

      // sim.mu_2            = 3*std::sqrt(3.0) / 3.0 * std::tan(32.76 * M_PI / 180);
      sim.mu_2            = std::tan(32.76 * M_PI / 180);

      //// For sintering only
      // sim.sinter_tc = 20;   // 1     // 20
      // sim.sinter_ec = 0.03; // 1e-3  // 0.03
      // sim.sinter_Sc = 20;
      // sim.particles.sinter_S.resize(sim.Np); std::fill( sim.particles.sinter_S.begin(), sim.particles.sinter_S.end(), sim.sinter_Sinf );

      // For MCC only:
      T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
      sim.particles.eps_pl_vol_mcc.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_mcc.begin(), sim.particles.eps_pl_vol_mcc.end(), eps_pl_vol_init );

      sim.simulate();
      // sim.validateRMA();


	return 0;
}
