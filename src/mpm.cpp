#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * 3D
//  * Cubic spline
//  * Complete doxygen
//  * Alembic output
//  * FLIP
//  * Stationary and moving levelsets with BC
//  * External body force dep on Lagrangian coordinates
//  * Laplace of splines and regularization
//  * Parallilize for loops over particles
//  * Loop only over grid points within 2dx of the particle position
//  * Make struct of std::vectors containging particle data, and the same for grid data
//  * Fix bug in the update of deformation gradient. See simple translation test case
//////////////////////////////////////////////////////////////////

int main(){
      std::cout << "This is Larsie" << std::endl;

      //csv2abc("test.abc");
      //return 0;

      Simulation sim;

      sim.end_frame = 20;
      sim.frame_dt = 1.0 / 1000.0;

      sim.gravity = TV2::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;

      sim.dx = 0.1;

      unsigned int Nloop = std::round(1.0/sim.dx);
      debug("Nloop           = ", Nloop);

      sim.Np = Nloop * Nloop * 4;

      sim.amplitude = 10.0;
      sim.neoHookean = true;
      sim.plasticity = false;
      sim.yield_stress = std::sqrt(2.0/3.0) * /* q_max */ 50000.0;

      sim.initialize(/* E */ 1e7, /* nu */ 0.3, /* rho */ 100);
      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass   = ", sim.particle_mass);
      debug("Np              = ", sim.Np);

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
                  T pvx = sim.amplitude;
                  T pvy = 0;
                  sim.particles_x(p) = px;
                  sim.particles_y(p) = py;
                  sim.particles_vx(p) = pvx;
                  sim.particles_vy(p) = pvy;
              } // end for d
          } // end for i
      } // end for j

      debug("Added particles = ", p);
      if ((p+1) != sim.Np){
          debug("Particle number mismatch!!!");
          return 0;
      }
    /////////////////////////////////////////////////////////////

    sim.simulate();

    /*

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
