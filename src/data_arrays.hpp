#ifndef DATAARRAYS_HPP
#define DATAARRAYS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include "tools.hpp"

class Particles{
public:
  Particles(unsigned int Np = 400){
      x  = TVX::Zero(Np);
      y  = TVX::Zero(Np);
      x0 = TVX::Zero(Np);
      y0 = TVX::Zero(Np);
      vx = TVX::Zero(Np);
      vy = TVX::Zero(Np);
      eps_pl_dev = TVX::Zero(Np);
      F.resize(Np); std::fill(F.begin(), F.end(), TM2::Identity());
  }
  TVX x;
  TVX y;
  TVX x0;
  TVX y0;
  TVX vx;
  TVX vy;
  TVX eps_pl_dev;
  std::vector<TM2> F;
};

class Grid{
public:
    Grid(unsigned int Nx = 21, unsigned int Ny = 21){
        x    = TVX::LinSpaced(Nx, -0.5, 1.5);
        y    = TVX::LinSpaced(Ny, -0.5, 1.5);
        vx   = TMX::Zero(Nx, Ny);
        vy   = TMX::Zero(Nx, Ny);
        mass = TMX::Zero(Nx, Ny);
    }
    TVX x;
    TVX y;
    TMX vx;
    TMX vy;
    TMX mass;
};


#endif  // DATAARRAYS_HPP
