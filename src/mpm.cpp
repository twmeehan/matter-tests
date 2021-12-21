#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * CUDA Parallilize
//////////////////////////////////////////////////////////////////

// REMEMBER TO SPECIFY DIMENSION IN TOOLS!!!!!

int main(){

      Simulation sim;
      sim.sim_name = "2d_thin";
      // sim.sim_name = "3d_elastic";
      // sim.sim_name = "test_rma_quad_anal";
      sim.directory = "/media/blatny/harddrive4/larsie/homomodel/";
      sim.end_frame = 800;
      sim.fps = 120;       // sim.fps = 8;
      sim.gravity = TV::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;
      sim.dt_max_coeff = 0.4;
      sim.flip_ratio = 0.99;
      sim.n_threads = 16;

      sim.initialize(/* E */ 1e8, /* nu */ 0.3, /* rho */ 300);

      ///////////// PARTICLES ON GRID //////////////////
      // sim.Lx = 1
      // sim.Ly = 2
      // int Npx = 30+1;
      // int Npy = 60+1;
      //
      // T dxp = sim.Lx / (Npx-1.0);
      // sim.dx = 2.0 * dxp;
      // sim.Np = Npx * Npy;// * Npz;
      ///////////////////////////////////////////////////

      ///////////// PARTICLES FROM FILE /////////////////
      // std::string sample = "samples/sample_Lx1.0_Ly2.0_r0.011/";
      std::string sample = "samples/sample_Lx20.0_Ly3.0_r0.05/";
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

      ///////////////////////////////////////////////////

      // for(int p = 0; p < sim.Np; p++){
      //     sim.particles.x[p][0] -= 1.0;
      // }

      sim.particle_volume = dx_p * dx_p;
      sim.particle_mass = sim.rho * sim.particle_volume;


      T vel_top = -0.2;       // T vel_top = -0.2/25.0;
      T vel_bot = 0.0;
      T vel_left = 0.0;
      T vel_right = 0.0;
     // T vel_back = 0.0;
     // T vel_front = 0.0;

      sim.vmin_factor = 25;
      sim.load_factor = 75;

      T offset = -0.01 * sim.dx/2.0; // When the grid is aligned with the boundary, it is important that the object overlap a bit into the particle domain
      T friction = 0.0;
      std::string name;
      // 3D
      // name = "Compressor"; InfinitePlate compressor = InfinitePlate(0, sim.Ly + offset, 0,    0, vel_top, 0, sim.vmin_factor, sim.load_factor,    top, SLIP, friction, name);  sim.objects.push_back(compressor);
      // name = "Ground";     InfinitePlate ground     = InfinitePlate(0, 0      - offset, 0,    0, vel_bot, 0,               1,               0, bottom, SLIP, friction, name);  sim.objects.push_back(ground);
      //
      // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0      - offset, 0, 0,    vel_left,  0, 0,             1,               0,   left, SLIP, friction, name);  sim.objects.push_back(side_left);
      // name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0, 0,    vel_right, 0, 0,             1,               0,  right, SLIP, friction, name);  sim.objects.push_back(side_right);
      //
      // name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0, 0, 0      - offset,    0, 0, vel_back,              1,               0,  back,  SLIP, friction, name);  sim.objects.push_back(side_back);
      // name = "SideFront";  InfinitePlate side_front = InfinitePlate(0, 0, sim.Lz + offset,    0, 0, vel_front,             1,               0,  front, SLIP, friction, name);  sim.objects.push_back(side_front);

      // 2D
      name = "Compressor"; InfinitePlate compressor = InfinitePlate(0, sim.Ly + offset,    0, vel_top, sim.vmin_factor, sim.load_factor,    top, SLIP, friction, name);  sim.objects.push_back(compressor);
      name = "Ground";     InfinitePlate ground     = InfinitePlate(0, 0      - offset,    0, vel_bot,               1,               0, bottom, SLIP, friction, name);  sim.objects.push_back(ground);

      name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0      - offset, 0,    vel_left,  0,             1,              0,   left, SLIP, friction, name);  sim.objects.push_back(side_left);
      name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx + offset, 0,    vel_right, 0,             1,              0,  right, SLIP, friction, name);  sim.objects.push_back(side_right);


      // Elastoplasticity
      sim.elastic_model = StvkWithHencky;
      // sim.plastic_model = VonMises;
      sim.plastic_model = Curved;
      // sim.plastic_model = NoPlasticity;
      sim.beta = 0.3; // 0.43-0.40*0.3;
      sim.M = 1.35;
      sim.p0 = 100e3; // 200e3;

      sim.xi = 1; //0.005;
      sim.xi_nonloc = 10; //0.5;

      sim.nonlocal_l = 0;

      T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
      sim.particles.eps_pl_vol_mcc.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_mcc.begin(), sim.particles.eps_pl_vol_mcc.end(), eps_pl_vol_init );


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
      // sim.validateRMA();\


	return 0;
}
