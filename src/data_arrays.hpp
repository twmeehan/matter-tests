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
      picx = TVX::Zero(Np);
      picy = TVX::Zero(Np);
      flipx = TVX::Zero(Np);
      flipy = TVX::Zero(Np);
      eps_pl_dev     = TVX::Zero(Np);
      eps_pl_vol     = TVX::Zero(Np);
      regularization = TVX::Zero(Np);
      cohesion_proj  = TVX::Zero(Np);
      tau.resize(Np); std::fill( tau.begin(), tau.end(), TM2::Zero()     );
      F.resize(Np);   std::fill( F.begin(),   F.end(),   TM2::Identity() );
  }
  TVX x;
  TVX y;
  TVX x0;
  TVX y0;
  TVX vx;
  TVX vy;
  TVX picx;
  TVX picy;
  TVX flipx;
  TVX flipy;
  TVX eps_pl_dev;
  TVX eps_pl_vol;
  TVX regularization;
  TVX cohesion_proj;
  std::vector<TM2> tau;
  std::vector<TM2> F;
};

class Grid{
public:
    Grid(unsigned int Nx = 21, unsigned int Ny = 21){
        x = TVX::LinSpaced(Nx, -0.5, 1.5);
        y = TVX::LinSpaced(Ny, -0.5, 1.5);
        mass    = TMX::Zero(Nx, Ny);
        vx      = TMX::Zero(Nx, Ny);
        vy      = TMX::Zero(Nx, Ny);
        flipx   = TMX::Zero(Nx, Ny);
        flipy   = TMX::Zero(Nx, Ny);
        regularization = TMX::Zero(Nx, Ny);
    }
    TVX x;
    TVX y;
    TMX vx;
    TMX vy;
    TMX flipx;
    TMX flipy;
    TMX mass;
    TMX regularization;
    T xc;
    T yc;
};


#endif  // DATAARRAYS_HPP
