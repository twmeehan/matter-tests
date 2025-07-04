#ifndef OBJECTBOXROTATED_HPP
#define OBJECTBOXROTATED_HPP

#include "object_general.hpp"

class ObjectBoxRotated : public ObjectGeneral {
public:
    TV center;
    TV half_extents;
    TM rotation;        // TM = matrix type, likely Eigen::Matrix<T, dim, dim>

    ObjectBoxRotated(BC bc_in, T friction_in, const TV& center_in, const TV& half_extents_in, const TM& rotation_in, std::string name_in = "")
        : ObjectGeneral(bc_in, friction_in, name_in),
          center(center_in),
          half_extents(half_extents_in),
          rotation(rotation_in) {}

    bool inside(const TV& X_world) const override {
        TV local = rotation.transpose() * (X_world - center); // World → local
        for (int i = 0; i < local.size(); ++i) {
            if (std::abs(local(i)) > half_extents(i)) return false;
        }
        return true;
    }

    TV normal(const TV& X_world) const override {
        TV local = rotation.transpose() * (X_world - center); // World → local

        int closest_axis = 0;
        T max_dist = 0.0;
        for (int i = 0; i < local.size(); ++i) {
            T dist = std::abs(half_extents(i) - std::abs(local(i)));
            if (dist > max_dist) {
                max_dist = dist;
                closest_axis = i;
            }
        }

        TV local_normal = TV::Zero();
        local_normal(closest_axis) = (local(closest_axis) > 0) ? 1 : -1;

        return rotation * local_normal; // Local → world
    }
};

#endif // OBJECTBOXROTATED_HPP
