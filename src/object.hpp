#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"


class GroundObject{
public:
  GroundObject(){
      vy_object = 0;
      vy_object = 0;
      y_object = -0.5;
  }

  bool inside(T x, T y){
      return distance(x, y) <= 0; // inside if dist is negative
  }

  T distance(T x, T y){
      return y - y_object;
  }

  void move(T dt){
      y_object += dt * vy_object;
  }

  T vx_object;
  T vy_object;
  T y_object;

};


#endif  // OBJECT_HPP
