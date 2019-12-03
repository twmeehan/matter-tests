#include "simulation.hpp"
#include "csv2abc.hpp"

// TODO:
// 	* Fix type template for double/float
//  * 3D
//  * Cubic spline
//  * adaptive timestepping
//  * complete doxygen

//////////////////////////////////////////////////////////////////

int main(){
    std::cout << "This is Larsie" << std::endl;

    csv2abc("test.abc");
    return 0;

    Simulation sim;

    //////////////////////////////////////////////////////////////
    ////////////// Unit box with 10 times dx = 0.1 ///////////////
    //////////////////////////////////////////////////////////////

      sim.Nt = 200;
      sim.dx = 0.1;
      sim.rho = 700;
      sim.setElasticParams(1e5, 0.3, sim.rho);
      debug("dt_max = ", sim.dt_max);

      sim.dt = 1e-6;
      debug("dt     = ", sim.dt);

      sim.Nx = 10;
      sim.Ny = 10;
      sim.Np = sim.Nx * sim.Ny * 4;
      debug("Np = ", sim.Np);

      sim.particle_volume = 1.0 / sim.Np;
      sim.particle_mass = sim.rho * 1.0 / sim.Np;

      sim.grid_X = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_Y = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_VX = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_VY = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);
      sim.grid_mass = Eigen::MatrixXd::Zero(sim.Nx, sim.Ny);

      Particle part;
      sim.particles.resize(sim.Np);
      std::fill(sim.particles.begin(), sim.particles.end(), part);

      sim.particles_x  = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_y  = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_vx = Eigen::VectorXd::Zero(sim.Np);
      sim.particles_vy = Eigen::VectorXd::Zero(sim.Np);

      double amplitude = 0.001;

      int p = -1;
      for(int i = 0; i < sim.Nx; i++){
          for(int j = 0; j < sim.Ny; j++){

              p++;
              double px = (i+0.25)*sim.dx;
              double py = (j+0.25)*sim.dx;
              double pvx = amplitude*std::sin(M_PI*px);
              double pvy = amplitude*std::sin(M_PI*py);
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.75)*sim.dx;
              py = (j+0.75)*sim.dx;
              pvx = amplitude*std::sin(M_PI*px);
              pvy = amplitude*std::sin(M_PI*py);
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.25)*sim.dx;
              py = (j+0.75)*sim.dx;
              pvx = amplitude*std::sin(M_PI*px);
              pvy = amplitude*std::sin(M_PI*py);
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;

              p++;
              px = (i+0.75)*sim.dx;
              py = (j+0.25)*sim.dx;
              pvx = amplitude*std::sin(M_PI*px);
              pvy = amplitude*std::sin(M_PI*py);
              sim.particles_x(p) = px;
              sim.particles_y(p) = py;
              sim.particles_vx(p) = pvx;
              sim.particles_vy(p) = pvy;
          } // end for i
      } // end for j

    /////////////////////////////////////////////////////////////

    sim.simulate();
	return 0;
}
