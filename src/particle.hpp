#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

typedef Eigen::Vector2d TV2;
typedef Eigen::Matrix2d TM2;

class Particle{
public:
  Particle(){
      F = TM2::Identity();
  }
  TM2 F;
};


#endif  // PARTICLE_HPP
