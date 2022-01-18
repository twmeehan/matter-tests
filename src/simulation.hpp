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
#include "timer.hpp"

#define OMP //comment out if not using omp

class Simulation{
public:
  Simulation();

#ifdef THREEDIM
  const unsigned int dim = 3;
#else
  const unsigned int dim = 2;
#endif

  std::string sim_name;
  std::string directory;

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

#ifdef THREEDIM
  T Lz;
#endif

  T rho;
  TV gravity;

  T amplitude;

  // Particle data
  unsigned int Np;
  T particle_mass;
  T particle_volume; // initial particle volume V0
  Particles particles;

  // Grid data
  unsigned int Nx, Ny;
#ifdef THREEDIM
  unsigned int Nz;
#endif

  unsigned int grid_nodes;
  Grid grid;

#ifdef THREEDIM
    inline unsigned int ind(unsigned int i, unsigned int j, unsigned int k){
      return (i*Ny + j) * Nz + k; // 3D
    }
#else
    inline unsigned int ind(unsigned int i, unsigned int j){
        return (i*Ny + j); // 3D
    }
#endif

  unsigned int vmin_factor;
  unsigned int load_factor;

  // Remeshing Fixed Grid
  T min_x_init;
  T max_x_init;
  T Nx_init;
  T low_x_init;
  T high_x_init;

  T min_y_init;
  T max_y_init;
  T Ny_init;
  T low_y_init;
  T high_y_init;

#ifdef THREEDIM
  T min_z_init;
  T max_z_init;
  T Nz_init;
  T low_z_init;
  T high_z_init;
#endif

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
  T yield_stress_min;

  // Drucker Prager
  T dp_slope;
  T dp_cohesion;

  // Perzyna
  T perzyna_exp;
  T perzyna_visc;

 // QuadraticLars
 T M;
 T beta;
 T p0;

  // Regularization by Laplacian
  T nonlocal_l;
  T nonlocal_l_sq;
  unsigned int nonlocal_support;

  // Objects
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
  void initialize(T E, T nu, T density);
  void simulate();
  void saveInfo();
  void saveParticleData(std::string extra = "");
  void saveGridData(std::string extra = "");

  void createDirectory();

  void advanceStep();
  void updateDt();
  void remesh();
  void remeshFixedInit();
  void remeshFixedCont();

  void P2G();
  void explicitEulerUpdate();
  void G2P();
  void deformationUpdate();

#ifdef OMP
      void P2G_Optimized_Parallel();
      void G2P_Optimized_Parallel();
      void deformationUpdate_Parallel();
      void explicitEulerUpdate_Optimized_Parallel();
      void plasticity_projection();
      void G2P_nonlocal();
      void P2G_nonlocal();
#else
    void P2G_Baseline();
    void P2G_Optimized();

    void explicitEulerUpdate_Baseline();
    void explicitEulerUpdate_Optimized();

    void G2P_Baseline();
    void G2P_Optimized();

    void deformationUpdate_Baseline();
#endif

  void positionUpdate();

  TM NeoHookeanPiola(TM & Fe);
  TM StvkWithHenckyPiola(TM & Fe);
  void plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial);

  void moveObjects();

#ifdef THREEDIM
  void boundaryCollision(T xi, T yi, T zi, TV& vi);
  void overwriteGridVelocity(T xi, T yi, T zi, TV& vi);
#else
  void boundaryCollision(T xi, T yi, TV& vi);
  void overwriteGridVelocity(T xi, T yi, TV& vi);
#endif

  void calculateMomentumOnParticles();
  void calculateMomentumOnGrid();
  void calculateMassConservation();

  void validateRMA();

  void periodicBoundaryConditions();
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
