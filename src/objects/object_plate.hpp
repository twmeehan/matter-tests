#ifndef OBJECTPLATE_HPP
#define OBJECTPLATE_HPP

#include "../tools.hpp"

class ObjectPlate{
public:

#ifdef THREEDIM

    ObjectPlate(T pos_object, T pos_upper, T pos_lower, PlateType plate_type, BoundaryCondition bc, T friction, std::string name, T vx_object, T vy_object, T vz_object, T vmin_factor, T load_factor) :
              pos_object(pos_object),
              pos_upper(pos_upper),
              pos_lower(pos_lower),
              plate_type(plate_type),
              bc(bc),
              friction(friction),
              name(name),
              vx_object(vx_object),
              vy_object(vy_object),
              vz_object(vz_object),
              vx_object_original(vx_object),
              vy_object_original(vy_object),
              vz_object_original(vz_object),
              vmin_factor(vmin_factor),
              load_factor(load_factor) {}

      bool inside(const TV& X_in){
          if (plate_type == left){
              if (X_in(1) < pos_upper && X_in(1) > pos_lower && (X_in(0) - pos_object) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          if (plate_type == right){
              if (X_in(1) < pos_upper && X_in(1) > pos_lower && (pos_object - X_in(0)) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          if (plate_type == bottom){
              if (X_in(0) < pos_upper && X_in(0) > pos_lower && (X_in(1) - pos_object) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          if (plate_type == top){
              if (X_in(0) < pos_upper && X_in(0) > pos_lower && (pos_object - X_in(1)) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          if (plate_type == back){
              if (X_in(2) - pos_object <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          if (plate_type == front){
              if (pos_object - X_in(2) <= 0){
                  return true;
              } else{
                  return false;
              }
          }

      } // end inside(x,y,z)


#else // TWODIM

    ObjectPlate(T pos_object, T pos_upper, T pos_lower, PlateType plate_type, BoundaryCondition bc, T friction, std::string name, T vx_object, T vy_object, T vmin_factor, T load_factor) :
              pos_object(pos_object),
              pos_upper(pos_upper),
              pos_lower(pos_lower),
              plate_type(plate_type),
              bc(bc),
              friction(friction),
              name(name),
              vx_object(vx_object),
              vy_object(vy_object),
              vx_object_original(vx_object),
              vy_object_original(vy_object),
              vmin_factor(vmin_factor),
              load_factor(load_factor) {}

    bool inside(const TV& X_in){
        if (plate_type == left){
            if (X_in(1) < pos_upper && X_in(1) > pos_lower && (X_in(0) - pos_object) <= 0){
                return true;
            } else{
                return false;
            }
        }
        if (plate_type == right){
            if (X_in(1) < pos_upper && X_in(1) > pos_lower && (pos_object - X_in(0)) <= 0){
                return true;
            } else{
                return false;
            }
        }
        if (plate_type == bottom){
            if (X_in(0) < pos_upper && X_in(0) > pos_lower && (X_in(1) - pos_object) <= 0){
                return true;
            } else{
                return false;
            }
        }
        if (plate_type == top){
            if (X_in(0) < pos_upper && X_in(0) > pos_lower && (pos_object - X_in(1)) <= 0){
                return true;
            } else{
                return false;
            }
        }
    }

#endif


  void move(T dt, T frame_dt, T time){

      T load_time = load_factor * frame_dt;
      if (time < load_time) {
          vx_object = vx_object_original / vmin_factor;
          vy_object = vy_object_original / vmin_factor;
          #ifdef THREEDIM
          vz_object = vz_object_original / vmin_factor;
          #endif
      }
      else{
          vx_object = vx_object_original;
          vy_object = vy_object_original;
          #ifdef THREEDIM
          vz_object = vz_object_original;
          #endif
      }

      if (plate_type == left || plate_type == right){
          pos_object += dt * vx_object;
          pos_upper  += dt * vy_object;
          pos_lower  += dt * vy_object;
      }
      else if (plate_type == bottom || plate_type == top){
          pos_object += dt * vy_object;
          pos_upper  += dt * vx_object;
          pos_lower  += dt * vx_object;
      }
#ifdef THREEDIM
      else if (plate_type == back || plate_type == front){
          pos_object += dt * vz_object;
      }
#endif


  } // end move(dt, frame_dt, time)

  T pos_object;
  T pos_upper;
  T pos_lower;

  BoundaryCondition bc;
  T friction;
  PlateType plate_type;
  std::string name;

  T vx_object;
  T vy_object;
  #ifdef THREEDIM
    T vz_object;
  #endif

  T vx_object_original;
  T vy_object_original;
  #ifdef THREEDIM
    T vz_object_original;
  #endif

  T vmin_factor;
  T load_factor;

};


#endif  // OBJECTPLATE_HPP
