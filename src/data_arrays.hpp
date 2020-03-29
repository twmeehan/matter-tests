#ifndef DATAARRAYS_HPP
#define DATAARRAYS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include "tools.hpp"

class Particles{
public:
  Particles(unsigned int Np = 400){
      x.resize(Np);    std::fill( x.begin(),    x.end(),    TV::Zero() );
      v.resize(Np);    std::fill( v.begin(),    v.end(),    TV::Zero() );
      pic.resize(Np);  std::fill( pic.begin(),  pic.end(),  TV::Zero() );
      flip.resize(Np); std::fill( flip.begin(), flip.end(), TV::Zero() );

      eps_pl_dev.resize(Np);     std::fill( eps_pl_dev.begin(),     eps_pl_dev.end(),     0.0 );
      eps_pl_vol.resize(Np);     std::fill( eps_pl_vol.begin(),     eps_pl_vol.end(),     0.0 );
      regularization.resize(Np); std::fill( regularization.begin(), regularization.end(), 0.0 );

      tau.resize(Np); std::fill( tau.begin(), tau.end(), TM::Zero()     );
      F.resize(Np);   std::fill( F.begin(),   F.end(),   TM::Identity() );
  }

  std::vector<TV> x;
  std::vector<TV> x0;
  std::vector<TV> v;
  std::vector<TV> pic;
  std::vector<TV> flip;

  std::vector<T> eps_pl_dev;
  std::vector<T> eps_pl_vol;
  std::vector<T> regularization;
  std::vector<T> cohesion_proj;

  std::vector<TM> tau;
  std::vector<TM> F;

};

class Grid{
public:
    Grid(unsigned int Nx = 21, unsigned int Ny = 21){
        x.resize(Nx); std::fill( x.begin(), x.end(), 0.0 );
        y.resize(Ny); std::fill( y.begin(), y.end(), 0.0 );

        v.resize(Nx*Ny);    std::fill( v.begin(),    v.end(),    TV::Zero() );
        flip.resize(Nx*Ny); std::fill( flip.begin(), flip.end(), TV::Zero() );

        mass.resize(Nx*Ny);           std::fill( mass.begin(),           mass.end(),           0.0 );
        regularization.resize(Nx*Ny); std::fill( regularization.begin(), regularization.end(), 0.0 );
    }
    std::vector<T> x;
    std::vector<T> y;
    std::vector<TV> v;
    std::vector<TV> flip;
    std::vector<T> mass;
    std::vector<T> regularization;
    T xc;
    T yc;
};


#endif  // DATAARRAYS_HPP
