#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "tools.hpp"
//#include "particle.hpp"

class Simulation{
public:
  //Simulation() : current_time_step(0), exit(0) {};
  Simulation();
  void initialize(T E, T nu, T density);
  void simulate();
  void saveSim(std::string extra = "");
  void saveGridVelocities(std::string extra = "");

  unsigned int current_time_step;
  unsigned int frame;
  unsigned int end_frame;
  T time;
  T final_time;
  T frame_dt;
  T dt;
  T dt_max;
  T wave_speed;
  T cfl;
  T dx;

  T rho;
  TV2 gravity;

  unsigned int Np;
  T particle_mass;
  T particle_volume; // initial particle volume V0

  std::vector<TM2> particles_F;
  TVX particles_x;
  TVX particles_y;
  TVX particles_x0;
  TVX particles_y0;
  TVX particles_vx;
  TVX particles_vy;

  unsigned int Nx, Ny;
  TVX lin_X;
  TVX lin_Y;
  TMX grid_VX;
  TMX grid_VY;
  TMX grid_mass;

  int exit;

  T mu;
  T lambda;

  T amplitude;

  bool neoHookean;
  bool plasticity;
  T yield_stress;

  void advanceStep();

  // advanceStep relies on (in order):
  void updateDt();
  void remesh();
  void P2G();
  void explicitEulerUpdate();
  void G2P();
  void deformationUpdate();
  void positionUpdate();

  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();
  std::pair<TMX, TMX> createExternalGravity();


}; // END Simulation class



#endif  // SIMULATION_HPP
