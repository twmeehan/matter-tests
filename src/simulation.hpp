#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include <vector>
#include <chrono>
#include "tools.hpp"
#include "object.hpp"
#include "data_arrays.hpp"
#include "timer.hpp"

#define OMP //comment out if not using omp

class Simulation{
public:
  //Simulation() : current_time_step(0), exit(0) {};
  Simulation();
  void initialize(T E, T nu, T density);
  void simulate();
  void saveParticleData(std::string extra = "");
  void saveGridData(std::string extra = "");

  std::string sim_name;

  unsigned int n_threads;
  unsigned int current_time_step;
  unsigned int frame;
  unsigned int end_frame;
  int exit;
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

  T amplitude;

  // Particle data
  unsigned int Np;
  T particle_mass;
  T particle_volume; // initial particle volume V0
  Particles particles;

  // Grid data
  unsigned int Nx, Ny;
  Grid grid;

  // Elastoplasticity
  T mu;
  T lambda;
  ElasticModel elastic_model;
  PlasticModel plastic_model;
  T yield_stress;
  T reg_length;
  T reg_const;
  // Objects
  BoundaryCondition boundary_condition;
  std::vector<InfinitePlate> objects;

  // Runtime measurements
  T runtime_p2g;
  T runtime_g2p;
  T runtime_euler;
  T runtime_defgrad;

  // Functions
  void advanceStep();
  void updateDt();
  void remesh();

  void P2G();
  void P2G_Baseline();
  void P2G_Optimized();

  void explicitEulerUpdate();
  void explicitEulerUpdate_Baseline();
  void explicitEulerUpdate_Optimized();

  void G2P();
  void G2P_Baseline();
  void G2P_Optimized();

  void deformationUpdate();
  void deformationUpdate_Baseline();

  #ifdef OMP
      void P2G_Optimized_Parallel();
      void G2P_Optimized_Parallel();
      void deformationUpdate_Parallel();
      void explicitEulerUpdate_Optimized_Parallel();
  #endif

  void positionUpdate();

  void moveObjects(T delta_t);
  void boundaryCollision(T xi, T yi, T& vxi, T& vyi);
  void boundaryCorrection(T xi, T yi, T& vxi, T& vyi);

  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();
  void calculateMassConservation();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();


}; // END Simulation class



#endif  // SIMULATION_HPP
