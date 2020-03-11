#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include "tools.hpp"
#include "object.hpp"
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

  // Particle data
  std::vector<TM2> particles_F;
  TVX particles_x;
  TVX particles_y;
  TVX particles_vx;
  TVX particles_vy;
  TVX particles_x0; // Lagrangian coord
  TVX particles_y0; // Lagrangian coord

  // Grid data
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
  unsigned int bc_type;

  bool neoHookean;
  bool plasticity;
  T yield_stress;

  GroundObject ground_object;

  void advanceStep();

  // advanceStep relies on (in order):
  void updateDt();
  void remesh();
  void P2G();
  void explicitEulerUpdate();
  void G2P();
  void deformationUpdate();
  void positionUpdate();

  void boundaryCollision(T xi, T yi, T& vxi, T& vyi);

  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();


}; // END Simulation class



#endif  // SIMULATION_HPP
