#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "tools.hpp"
#include "object.hpp"
#include "data_arrays.hpp"
#include "timer.hpp"

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
  Particles particles;

  // Grid data
  unsigned int Nx, Ny;
  Grid grid;

  int exit;

  T mu;
  T lambda;

  T amplitude;
  unsigned int bc_type;

  bool neoHookean;
  bool plasticity;
  T yield_stress;

  std::vector<InfinitePlate> objects;

  T runtime_p2g;
  T runtime_g2p;
  T runtime_euler;
  T runtime_defgrad;

  void advanceStep();
  void updateDt();
  void remesh();
  void P2G();
  void P2G_Optimized();
  void P2G_Baseline();
  void explicitEulerUpdate();
  void G2P();
  void G2P_Optimized();
  void G2P_Baseline();
  void deformationUpdate();
  void positionUpdate();

  void moveObjects();
  void boundaryCollision(T xi, T yi, T& vxi, T& vyi);

  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();
  void calculateMassConservation();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();


}; // END Simulation class



#endif  // SIMULATION_HPP
