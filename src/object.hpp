#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"
#include <string>

class InfinitePlate{
public:
  InfinitePlate(T x_object, T y_object, T z_object, T vx_object, T vy_object, T vz_object, T vmin_factor, T load_factor, PlateType plate_type, BoundaryCondition bc, std::string name) :
                x_object(x_object),
                y_object(y_object),
                z_object(z_object),
                vx_object(vx_object),
                vy_object(vy_object),
                vz_object(vz_object),
                vx_object_original(vx_object),
                vy_object_original(vy_object),
                vz_object_original(vz_object),
                vmin_factor(vmin_factor),
                load_factor(load_factor),
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

      T load_time = load_factor * frame_dt;
      if (time < load_time) {
          vx_object = vx_object_original / vmin_factor;
          vy_object = vy_object_original / vmin_factor;
          vz_object = vz_object_original / vmin_factor;
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

  T vmin_factor;
  T load_factor;

  BoundaryCondition bc;
  PlateType plate_type;
  std::string name;

};


// class TopPlate : public InfinitePlate{
// public:
//   TopPlate(T y_object, T vy_object, BoundaryCondition bc) : InfinitePlate(0, y_object, 0, 0, vy_object, 0, top, bc, "TopPlate") {}
//
//   bool inside(T x, T y, T z){
//       return distance(x,y,z) <= 0; // inside if dist is negative
//   }
//
//   T distance(T x, T y, T z) {
//       return (y_object - y);
//   }
//
//   void move(T dt, T frame_dt, T time) {
//
//       if (time < frame_dt) {
//           vy_object = (vy_object_original / frame_dt) * time;
//       }
//       else{
//           vy_object = vy_object_original;
//       }
//       y_object += dt * vy_object;
//   }
//
// };


#endif  // OBJECT_HPP
