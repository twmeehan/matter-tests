#include "simulation.hpp"
//#include "csv2abc.hpp"

// TODO:
// 	* Fix type template for double/float
//  * 3D
//  * Cubic spline
//  * adaptive timestepping
//  * complete doxygen

//////////////////////////////////////////////////////////////////

int main(){
    std::cout << "This is Larsie" << std::endl;

    //csv2abc("test.abc");
    //return 0;

    Simulation sim;

    //////////////////////////////////////////////////////////////
    ////////////// Unit box with 10 times dx = 0.1 ///////////////
    //////////////////////////////////////////////////////////////

      sim.T = 0.02;
      sim.max_time_steps = 500;
      sim.dx = 0.1;
      sim.rho = 1000;
      sim.setElasticParams(1e7, 0.3, sim.rho);
      debug("dt_max = ", sim.dt_max);

      sim.cfl = 0.6;
      sim.dt = 1e-3;
      debug("dt     = ", sim.dt);

      // N = (L / dx + 1)
      sim.Nx = 21;
      sim.Ny = 21;

      ////////////// If remesh every time step /////////////////
      // sim.grid_X = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      // sim.grid_Y = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      // sim.grid_VX = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      // sim.grid_VY = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      // sim.grid_mass = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);

      ///////////////// If constant mesh ////////////////////////
      Eigen::VectorXd lin_x = Eigen::VectorXd::LinSpaced(sim.Nx, -0.5, 1.5);
      Eigen::VectorXd lin_y = Eigen::VectorXd::LinSpaced(sim.Ny, -0.5, 1.5);

      sim.grid_X.resize(sim.Ny, sim.Nx);
      sim.grid_Y.resize(sim.Ny, sim.Nx);
      for (int i = 0; i < sim.Ny; ++i) {
        sim.grid_X.row(i) = lin_x.transpose();
      }
      for (int j = 0; j < sim.Nx; ++j) {
        sim.grid_Y.col(j) = lin_y;
      }
      sim.grid_X.transposeInPlace();
      sim.grid_Y.transposeInPlace();

      sim.grid_VX   = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_VY   = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_mass = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      ////////////////////////////////////////////////////////////////

      sim.Np = 10 * 10 * 4;
      debug("Np = ", sim.Np);

      sim.particle_volume = 1.0 / sim.Np; // INITIAL particle volume V^0
      sim.particle_mass = sim.rho * sim.particle_volume;

      debug("particle_volume = ", sim.particle_volume);
      debug("particle_mass = ", sim.particle_mass);


      Particle part;
      sim.particles.resize(sim.Np);
      std::fill(sim.particles.begin(), sim.particles.end(), part);

      sim.particles_x  = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_y  = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_vx = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_vy = Eigen::VectorXd::Zero(sim.Np);

      double amplitude = 50;

      int p = -1;
      for(int i = 0; i < 10; i++){
          for(int j = 0; j < 10; j++){

              p++;
              double px = (i+0.25)*sim.dx;
              double py = (j+0.25)*sim.dx;
              double pvx = amplitude*std::sin( M_PI*(px-0.5) );
              double pvy = amplitude*std::sin( M_PI*(py-0.5) );
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

    /////////////////////////////////////////////////////////////

//    sim.simulate();


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

/*
NB!!!!!!!!
THE MISTAKE IS EITHER IN explicitEulerUpdate() OR deformationUpdate()
*/

	return 0;
}
