// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::boundaryCollision(int index, TV Xi, TV& vi){

    // Make a copy
    TV vi_orig = vi;

    // New positions
    Xi += dt * vi; // recall Xi was passed by value

    for(auto obj : objects) {
        bool colliding = obj->inside(Xi);
        if (colliding) {
            TV v_rel = vi_orig;

            if (obj->bc == STICKY) {
                v_rel.setZero();
            } // end STICKY

            else if (obj->bc == SEPARATE) {
                TV n = obj->normal(Xi);
                T dot = v_rel.dot(n);
                if (dot < 0){ // if moving towards object

                    T friction = obj->friction;
                    if (use_material_fricton)
                        friction = grid.friction[index];

                    TV v_tang = v_rel - dot * n;
                    if (friction > 0){
                        if( -dot * friction < v_tang.norm() )
                            v_rel = v_tang + v_tang.normalized() * dot * friction;
                        else
                            v_rel.setZero();
                    } else{
                        v_rel = v_tang;
                    } // end non-zero friction
                }

            } // end SEPARATE
            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi = v_rel;
        } // end if colliding

    } // end iterator over general objects


#ifdef THREEDIM

    for (ObjectPlate &obj : plates) {
        bool colliding = obj.inside(Xi);
        if (colliding) {
            T vx_rel = vi_orig(0) - obj.vx_object;
            T vy_rel = vi_orig(1) - obj.vy_object;
            T vz_rel = vi_orig(2) - obj.vz_object;

            if (obj.bc == STICKY) {
                vx_rel = 0;
                vy_rel = 0;
                vz_rel = 0;
            } // end STICKY

            else if (obj.bc == SLIP) {
                if (obj.plate_type == top || obj.plate_type == bottom){
                    // tangential velocity is the (x,z) components

                    T vel_t      = std::sqrt(vx_rel*vx_rel + vz_rel*vz_rel);
                    T fric_vel_n = obj.friction * std::abs(vy_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vx_rel*vx_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vz_rel*vz_rel));
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end top and bottom plate
                else if (obj.plate_type == left || obj.plate_type == right){
                    // tangential velocity is the (y,z) components

                    T vel_t      = std::sqrt(vy_rel*vy_rel + vz_rel*vz_rel);
                    T fric_vel_n = obj.friction * std::abs(vx_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vy_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vy_rel*vy_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vz_rel*vz_rel));
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end left or right plate
                else if (obj.plate_type == front || obj.plate_type == back){
                    // tangential velocity is the (x,y) components

                    T vel_t      = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel);
                    T fric_vel_n = obj.friction * std::abs(vz_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vy_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vx_rel*vx_rel));
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vy_rel*vy_rel));
                    }

                    // normal component (z) must be set to zero
                    vz_rel = 0;

                } // end front or back plate
            } // end SLIP



            else if (obj.bc == SEPARATE) {
                if ((obj.plate_type == top && vy_rel > 0) || (obj.plate_type == bottom && vy_rel < 0)){
                    // tangential velocity is the (x,z) components
                    T vel_t      = std::sqrt(vx_rel*vx_rel + vz_rel*vz_rel);
                    T fric_vel_n = obj.friction * std::abs(vy_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vx_rel*vx_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vz_rel*vz_rel));
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end top and bottom plate
                else if ((obj.plate_type == left && vx_rel < 0) || (obj.plate_type == right && vx_rel > 0)){
                    // tangential velocity is the (y,z) components
                    T vel_t      = std::sqrt(vy_rel*vy_rel + vz_rel*vz_rel);
                    T fric_vel_n = obj.friction * std::abs(vx_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vy_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vy_rel*vy_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vz_rel*vz_rel));
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end left or right plate
                else if ((obj.plate_type == front && vz_rel > 0) || (obj.plate_type == back && vz_rel < 0 )){
                    // tangential velocity is the (x,y) components
                    T vel_t      = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel);
                    T fric_vel_n = obj.friction * std::abs(vz_rel);

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vy_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vx_rel*vx_rel));
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vy_rel*vy_rel));
                    }

                    // normal component (z) must be set to zero
                    vz_rel = 0;

                } // end front or back plate
            } // end SEPARATE



            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi(0) = vx_rel + obj.vx_object;
            vi(1) = vy_rel + obj.vy_object;
            vi(2) = vz_rel + obj.vz_object;
        } // end if colliding

    } // end iterator over 3D plate objects


#else // TWODIM


    for (ObjectPlate &obj : plates) {
        bool colliding = obj.inside(Xi);
        if (colliding) {
            T vx_rel = vi_orig(0) - obj.vx_object;
            T vy_rel = vi_orig(1) - obj.vy_object;

            if (obj.bc == STICKY) {
                vx_rel = 0;
                vy_rel = 0;
            } // end STICKY

            else if (obj.bc == SLIP) {
                if (obj.plate_type == top || obj.plate_type == bottom){
                    // tangential velocity is the (x) component and must be changed
                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < obj.friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= obj.friction * std::abs(vy_rel);
                        else
                            vx_rel += obj.friction * std::abs(vy_rel);
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end top and bottom plate
                else if (obj.plate_type == left || obj.plate_type == right){
                    // tangential velocity is the (y) component and must be changed

                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < obj.friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= obj.friction * std::abs(vx_rel);
                        else
                            vy_rel += obj.friction * std::abs(vx_rel);
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end left or right plate

            } // end SLIP



            else if (obj.bc == SEPARATE) {

                if ((obj.plate_type == top && vy_rel > 0) || (obj.plate_type == bottom && vy_rel < 0)){
                    // tangential velocity is the (x) component and must be changed
                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < obj.friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= obj.friction * std::abs(vy_rel);
                        else
                            vx_rel += obj.friction * std::abs(vy_rel);
                    }

                    // normal component (y) must be set to zero IF MOVING TOWARDS THE PLATE
                    vy_rel = 0;

                } // end top and bottom plate
                else if ((obj.plate_type == left && vx_rel < 0) || (obj.plate_type == right && vx_rel > 0)){
                    // tangential velocity is the (y) component and must be changed
                    if (obj.friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < obj.friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= obj.friction * std::abs(vx_rel);
                        else
                            vy_rel += obj.friction * std::abs(vx_rel);
                    }

                    // normal component (x) must be set to zero IF MOVING TOWARDS THE PLATE
                    vx_rel = 0;

                } // end left or right plate
            } // end SEPARATE



            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi(0) = vx_rel + obj.vx_object;
            vi(1) = vy_rel + obj.vy_object;
        } // end if colliding

    } // end iterator over 2D plate objects

#endif

} // end boundaryCollision
