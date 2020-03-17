#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * 3D
//  * Complete doxygen
//  * Alembic output
//  * FLIP
//  * Boundary condition for Laplace? Basically using Dirichlet now
//  * Parallilize explicit_euler_update
//////////////////////////////////////////////////////////////////

int main(){

      Simulation sim;

      sim.sim_name = "elastic_slip";

      sim.end_frame = 40;
      sim.frame_dt = 1.0 / 400.0;
      sim.dx = 0.05;
      sim.gravity = TV2::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;

      sim.n_threads = 4;

      unsigned int Nloop = std::round(1.0/sim.dx);
      debug("Nloop           = ", Nloop);
      sim.Np = Nloop * Nloop * 4;

      std::string name;
      name = "Ground";     InfinitePlate ground      = InfinitePlate(0, 0, 0, 0,  bottom, name); sim.objects.push_back(ground);
      name = "Compressor"; InfinitePlate compressor  = InfinitePlate(0, 1, 0, -1, top,    name); sim.objects.push_back(compressor);
      // name = "Left";       InfinitePlate left_plate  = InfinitePlate(0, 0, 0, 0,  left,   name); sim.objects.push_back(left_plate);
      // name = "Right";      InfinitePlate right_plate = InfinitePlate(1, 0, 0, 0,  right,  name); sim.objects.push_back(right_plate);

      sim.boundary_condition = STICKY;
      sim.friction = 0.5;

      sim.initialize(/* E */ 1e7, /* nu */ 0.0, /* rho */ 100);
      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass   = ", sim.particle_mass);
      debug("Np              = ", sim.Np);

      sim.elastic_model = StvkWithHencky;
      sim.plastic_model = DPSimpleSoft;

      // Von Mises:
      sim.yield_stress = std::sqrt(2.0/3.0) * /* q_max */ 400000.0;

      // DPSimpleSoft
      sim.friction_angle = 13;
      sim.cohesion = 10000.0 / (sim.K * sim.dim); // p_min = - K * dim * cohesion
      sim.xi = 1e15;

      // Regularization by Laplacian
      sim.reg_length = 0.00;
      sim.reg_const = 2*sim.mu;

      sim.amplitude = 0.0;
      std::vector<T> disp_i(4); disp_i[0] = 0.25; disp_i[1] = 0.75; disp_i[2] = 0.25; disp_i[3] = 0.75;
      std::vector<T> disp_j(4); disp_j[0] = 0.25; disp_j[1] = 0.75; disp_j[2] = 0.75; disp_j[3] = 0.25;
      int p = -1;
      for(int i = 0; i < Nloop; i++){
          for(int j = 0; j < Nloop; j++){
              for(int d = 0; d < 4; d++){
                  p++;
                  T px = (i+disp_i[d])*sim.dx;
                  T py = (j+disp_j[d])*sim.dx;
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
                  T pvx = 0;
                  T pvy = sim.amplitude;
                  sim.particles.x(p) = px;
                  sim.particles.y(p) = py;
                  sim.particles.vx(p) = pvx;
                  sim.particles.vy(p) = pvy;
              } // end for d
          } // end for i
      } // end for j

      debug("Added particles = ", p);
      if ((p+1) != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }


    sim.simulate();

    /////////////////////////////////////////////////////////////

    /*

    ///////////// ALEMBIC TESTING: ////////////////
    // csv2abc("test.abc");
    // return 0;

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
        // for (int i=0; i<50; i++){
        //     sim.updateDt();
        //     sim.P2G();
        //     sim.saveGridVelocities("before_");
        //     sim.explicitEulerUpdate();
        //     sim.saveGridVelocities("after_");
        //     sim.G2P();
        //     sim.saveSim("before_");
        //     sim.deformationUpdate();
        //     sim.positionUpdate();
        //     sim.saveSim("after_");
        //     sim.current_time_step++;
        // }



	return 0;
}
