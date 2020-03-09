#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
//  * 3D
//  * Cubic spline
//  * Complete doxygen
//  * Alembic output
//  * Plastic model
//  * Stationary and moving levelsets with BC
//  * Fix St Venant Kirchhoff with Hencky strain

//////////////////////////////////////////////////////////////////

int main(){
      std::cout << "This is Larsie" << std::endl;

      //csv2abc("test.abc");
      //return 0;

      Simulation sim;

      sim.end_frame = 20;
      sim.frame_dt = 0.0001;

      sim.gravity = TV2::Zero(); sim.gravity[1] = 0;
      sim.cfl = 0.6;

      sim.dx = 0.1;
      sim.Np = 10 * 10 * 4;

      sim.initialize(/* E */ 1e7, /* nu */ 0.3, /* rho */ 100);
      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass   = ", sim.particle_mass);
      debug("Np              = ", sim.Np);

      T amplitude = 100.0;

      int p = -1;
      for(int i = 0; i < 10; i++){
          for(int j = 0; j < 10; j++){

              p++;
              T px = (i+0.25)*sim.dx;
              T py = (j+0.25)*sim.dx;
              T pvx = amplitude*std::sin( M_PI*(px-0.5) );
              T pvy = amplitude*std::sin( M_PI*(py-0.5) );
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.75)*sim.dx;
              py = (j+0.75)*sim.dx;
              pvx = amplitude*std::sin( M_PI*(px-0.5) );
              pvy = amplitude*std::sin( M_PI*(py-0.5) );
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.25)*sim.dx;
              py = (j+0.75)*sim.dx;
              pvx = amplitude*std::sin( M_PI*(px-0.5) );
              pvy = amplitude*std::sin( M_PI*(py-0.5) );
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.75)*sim.dx;
              py = (j+0.25)*sim.dx;
              pvx = amplitude*std::sin( M_PI*(px-0.5) );
              pvy = amplitude*std::sin( M_PI*(py-0.5) );
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;
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
    NB!!!!!!!!
    THERE IS MISTAKE IN explicitEulerUpdate() for St Venant Kirchhoff with Hencky strain
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
