#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"
#include <string>

class InfinitePlate{
public:
  InfinitePlate() : x_object(0.0), y_object(0.0), vx_object(0.0), vy_object(0.0), plate_type(bottom), name("noname") {}
  InfinitePlate(T x_object, T y_object, T vx_object, T vy_object, PlateType plate_type, std::string name) : x_object(x_object), y_object(y_object), vx_object(vx_object), vy_object(vy_object), plate_type(plate_type), name(name) {}

  bool inside(T x, T y){
      return distance(x, y) <= 0; // inside if dist is negative
  }

  T distance(T x, T y){
      if (plate_type == bottom)
        return (y - y_object);
      else if (plate_type == top)
          return (y_object - y);
      else if (plate_type == left)
          return (x - x_object);
      else // right
          return (x_object - x);
  }

  void move(T dt){
      x_object += dt * vx_object;
      y_object += dt * vy_object;
  }

  T x_object;
  T y_object;

  T vx_object;
  T vy_object;

  PlateType plate_type;
  std::string name;

};


//
// class InfinitePlate{
// public:
//   InfinitePlate() : x_object(0.0), y_object(0.0), vx_object(0.0), vy_object(0.0), name("InfinitePlate") {}
//   InfinitePlate(std::string name) : x_object(0.0), y_object(0.0), vx_object(0.0), vy_object(0.0), name(name) {}
//   InfinitePlate(T x_object, T y_object, T vx_object, T vy_object, std::string name) : x_object(x_object), y_object(y_object), vx_object(vx_object), vy_object(vy_object), name(name) {}
//
//   bool inside(T x, T y){
//         return distance(x, y) <= 0; // inside if dist is negative
//   }
//
//   virtual T distance(T x, T y){ return nan(""); }
//
//   void move(T dt){
//       x_object += dt * vx_object;
//       y_object += dt * vy_object;
//   }
//
//   T x_object;
//   T y_object;
//   T vx_object;
//   T vy_object;
//   std::string name;
//
// };
//
//
//
//
// class TopPlate : public InfinitePlate{
// public:
//   TopPlate() : InfinitePlate("TopPlate") {}
//   TopPlate(T y_object, T vy_object) : InfinitePlate(0, y_object, 0, vy_object, "TopPlate") {}
//
//   T distance(T x, T y) override {
//       return (y_object - y);
//   }
//
// };
//
//
//
// class BottomPlate : public InfinitePlate{
// public:
//   BottomPlate() : InfinitePlate("BottomPlate") {}
//   BottomPlate(T y_object, T vy_object) : InfinitePlate(0, y_object, 0, vy_object, "BottomPlate") {}
//
//   T distance(T x, T y) override {
//       return (y - y_object);
//   }
//
// };


#endif  // OBJECT_HPP
