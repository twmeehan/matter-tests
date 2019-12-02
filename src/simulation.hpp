#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "tools.hpp"
#include "particle.hpp"

typedef Eigen::Vector2d TV2;
typedef Eigen::Matrix2d TM2;

class Simulation{
public:
  Simulation() : current_step(0), exit(0) {};
  void setElasticParams(double E, double nu, double density);
  void simulate();
  void saveSim();
  void saveGridVelocities();

  unsigned int Nt;
  unsigned int current_step;
  double dt;
  double dt_max;
  double dx;

  double rho;

  unsigned int Np;
  double particle_mass;
  double particle_volume; // initial particle volume V0

  std::vector<Particle> particles;
  Eigen::VectorXd particles_x;
  Eigen::VectorXd particles_y;
  Eigen::VectorXd particles_vx;
  Eigen::VectorXd particles_vy;

  unsigned int Nx, Ny;
  Eigen::MatrixXd grid_X;
  Eigen::MatrixXd grid_Y;
  Eigen::MatrixXd grid_VX;
  Eigen::MatrixXd grid_VY;
  Eigen::MatrixXd grid_mass;

private:

  int exit;

  double mu;
  double lambda;

  void advanceStep();

  // advanceStep relies on (in order):
  void remesh();
  void P2G();
  void explicitEulerUpdate();
  void G2P();
  void deformationUpdate();
  void positionUpdate();

}; // END Simulation class



#endif  // SIMULATION_HPP
