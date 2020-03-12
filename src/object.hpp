#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"
#include <string>

enum PlateType { upper, lower};

class InfinitePlate{
public:
  InfinitePlate() : y_object(0.0), vy_object(0.0), vx_object(0.0), plate_type(lower), name("noname") {}
  InfinitePlate(T y_object, T vy_object, PlateType plate_type, std::string name) : y_object(y_object), vy_object(vy_object), vx_object(0.0), plate_type(plate_type), name(name) {}

  bool inside(T x, T y){
      return distance(x, y) <= 0; // inside if dist is negative
  }

  T distance(T x, T y){
      if (plate_type == lower){
        return (y - y_object);
      }
      return (y_object - y);
  }

  void move(T dt){
      y_object += dt * vy_object;
  }

  T vx_object;
  T vy_object;
  T y_object;
  PlateType plate_type;
  std::string name;

};


#endif  // OBJECT_HPP
