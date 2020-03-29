#ifndef DATAARRAYS_HPP
#define DATAARRAYS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include "tools.hpp"

class Particles{
public:
  Particles(unsigned int Np = 400){
      x.resize(Np);    std::fill( x.begin(),    x.end(),    TV2::Zero() );
      //x0.resize(Np);   std::fill( x0.begin(),   x0.end(),   TV2::Zero() );
      v.resize(Np);    std::fill( v.begin(),    v.end(),    TV2::Zero() );
      pic.resize(Np);  std::fill( pic.begin(),  pic.end(),  TV2::Zero() );
      flip.resize(Np); std::fill( flip.begin(), flip.end(), TV2::Zero() );

      eps_pl_dev.resize(Np);     std::fill( eps_pl_dev.begin(),     eps_pl_dev.end(),     0.0 );
      eps_pl_vol.resize(Np);     std::fill( eps_pl_vol.begin(),     eps_pl_vol.end(),     0.0 );
      regularization.resize(Np); std::fill( regularization.begin(), regularization.end(), 0.0 );
      //cohesion_proj.resize(Np);  std::fill( cohesion_proj.begin(),  cohesion_proj.end(),  0.0 );

      tau.resize(Np); std::fill( tau.begin(), tau.end(), TM2::Zero()     );
      F.resize(Np);   std::fill( F.begin(),   F.end(),   TM2::Identity() );
  }

  std::vector<TV2> x;
  std::vector<TV2> x0;
  std::vector<TV2> v;
  std::vector<TV2> pic;
  std::vector<TV2> flip;

  std::vector<T> eps_pl_dev;
  std::vector<T> eps_pl_vol;
  std::vector<T> regularization;
  std::vector<T> cohesion_proj;

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
