#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"
#include <string>

class InfinitePlate{
public:
  InfinitePlate() :
                x_object(0.0),
                y_object(0.0),
                z_object(0.0),
                vx_object(0.0),
                vy_object(0.0),
                vz_object(0.0),
                vx_object_original(0.0),
                vy_object_original(0.0),
                vz_object_original(0.0),
                plate_type(bottom),
                bc(STICKY),
                name("noname") {}
  InfinitePlate(T x_object, T y_object, T z_object, T vx_object, T vy_object, T vz_object, PlateType plate_type, BoundaryCondition bc, std::string name) :
                x_object(x_object),
                y_object(y_object),
                z_object(z_object),
                vx_object(vx_object),
                vy_object(vy_object),
                vz_object(vz_object),
                vx_object_original(vx_object),
                vy_object_original(vy_object),
                vz_object_original(vz_object),
                plate_type(plate_type),
                bc(bc),
                name(name) {}

  bool inside(T x, T y, T z){
      return distance(x, y, z) <= 0; // inside if dist is negative
  }

  T distance(T x, T y, T z){
      if (plate_type == bottom)
        return (y - y_object);
      else if (plate_type == top)
          return (y_object - y);
      else if (plate_type == left)
          return (x - x_object);
      else if (plate_type == right)
          return (x_object - x);
      else if (plate_type == back)
          return (z - z_object);
      else if (plate_type == front)
          return (z_object - z);
  }

  void move(T dt, T frame_dt, T time){

      if (time < frame_dt) {
          vx_object = (vx_object_original / frame_dt) * time;
          vy_object = (vy_object_original / frame_dt) * time;
          vz_object = (vz_object_original / frame_dt) * time;
      }
      else{
          vx_object = vx_object_original;
          vy_object = vy_object_original;
          vz_object = vz_object_original;
      }

      x_object += dt * vx_object;
      y_object += dt * vy_object;
      z_object += dt * vz_object;
  }

  T x_object;
  T y_object;
  T z_object;

  T vx_object;
  T vy_object;
  T vz_object;

  T vx_object_original;
  T vy_object_original;
  T vz_object_original;

  BoundaryCondition bc;
  PlateType plate_type;
  std::string name;

};

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
