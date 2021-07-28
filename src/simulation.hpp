#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <assert.h>

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
  std::string directory;

  unsigned int dim;
  unsigned int n_threads;
  unsigned int current_time_step;
  unsigned int frame;
  unsigned int end_frame;
  int exit;
  T time;
  T final_time;
  T fps;
  T frame_dt;
  T dt;
  T dt_max;
  T dt_max_coeff;
  T flip_ratio;
  T wave_speed;
  T cfl;
  T dx;
  T Lx;
  T Ly;
  T Lz;
  T rho;
  TV gravity;

  T amplitude;

  // Particle data
  unsigned int Np;
  T particle_mass;
  T particle_volume; // initial particle volume V0
  Particles particles;

  // Grid data
  unsigned int Nx, Ny, Nz;
  Grid grid;
  inline unsigned int ind(unsigned int i, unsigned int j, unsigned int k){
      return (i*Ny + j) * Nz + k; // 3D
  }

  // Remeshing Fixed Grid
  T max_y_init;
  T Ny_init;
  T low_y_init;
  T high_y_init;

  // Elastoplasticity
  T mu;
  T lambda;
  T K;
  ElasticModel elastic_model;
  PlasticModel plastic_model;
  T xi;
  T xi_nonloc;

  // Von Mises:
  T yield_stress_orig;

  // DPSimpleSoft
  T friction_angle;
  T alpha_K_d_over_2mu;
  T cohesion;

 // QuadraticLars
 T M;
 T beta;
 T p0;

  // Regularization by Laplacian
  T nonlocal_l;
  T nonlocal_l_sq;
  unsigned int nonlocal_support;

  // Objects
  T friction;
  std::vector<InfinitePlate> objects;

  // Runtime measurements
  T runtime_p2g;
  T runtime_g2p;
  T runtime_euler;
  T runtime_defgrad;

  // Precomputations
  T one_over_dx;
  T one_over_dx_square;

  T mu_sqrt6;

  // Functions
  void createDirectory();

  void advanceStep();
  void updateDt();
  void remesh();
  void remeshFixedInit();
  void remeshFixedCont();

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
      void plasticity_projection();
      void G2P_nonlocal();
      void P2G_nonlocal();
  #endif

  void positionUpdate();

  TM NeoHookeanPiola(TM & Fe);
  TM StvkWithHenckyPiola(TM & Fe);
  void plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial);

  void moveObjects();
  void boundaryCollision(T xi, T yi, T zi, TV& vi);
  // void boundaryCorrection(T xi, T yi, T& vxi, T& vyi);
  void overwriteGridVelocity(T xi, T yi, T zi, TV& vi);
  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();
  void calculateMassConservation();

  void validateRMA();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();


}; // END Simulation class



inline TM Simulation::NeoHookeanPiola(TM & Fe){
    return mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
} // end NeoHookeanPiola

inline TM Simulation::StvkWithHenckyPiola(TM & Fe){
    Eigen::JacobiSVD<TM> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
    TA sigma = svd.singularValues().array(); // abs() for inverse also??
    TM logSigma = sigma.abs().log().matrix().asDiagonal();
    TM invSigma = sigma.inverse().matrix().asDiagonal();
    TM dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
    return dPsidF;
} // end StvkWithHenckyPiola



#endif  // SIMULATION_HPP
