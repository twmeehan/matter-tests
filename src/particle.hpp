#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include "tools.hpp"

class Particle{
public:
  Particle(){
      F = TM2::Identity();
  }
  TM2 F;
};


#endif  // PARTICLE_HPP
