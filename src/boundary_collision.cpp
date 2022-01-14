#include "simulation.hpp"

#ifdef THREEDIM

void Simulation::boundaryCollision(T xi, T yi, T zi, TV& vi){

    // Reference to velocity components
    T& vxi = vi(0);
    T& vyi = vi(1);
    T& vzi = vi(2);

    // New positions
    xi += dt * vxi;
    yi += dt * vyi;
    zi += dt * vzi;

    for (InfinitePlate &obj : objects) {
        bool colliding = obj.inside(xi, yi, zi);
        if (colliding) {
            T vx_rel = vxi - obj.vx_object;
            T vy_rel = vyi - obj.vy_object;
            T vz_rel = vzi - obj.vz_object;

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

            vxi = vx_rel + obj.vx_object;
            vyi = vy_rel + obj.vy_object;
            vzi = vz_rel + obj.vz_object;
        } // end if colliding

    } // end iterator over objects

} // end boundaryCollision

#else

void Simulation::boundaryCollision(T xi, T yi, TV& vi){

    // TODO: Make this a vector-based function

    // Reference to velocity components
    T& vxi = vi(0);
    T& vyi = vi(1);

    // New positions
    xi += dt * vxi;
    yi += dt * vyi;

    for (InfinitePlate &obj : objects) {
        bool colliding = obj.inside(xi, yi);
        if (colliding) {
            T vx_rel = vxi - obj.vx_object;
            T vy_rel = vyi - obj.vy_object;

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

            vxi = vx_rel + obj.vx_object;
            vyi = vy_rel + obj.vy_object;
        } // end if colliding

    } // end iterator over objects

} // end boundaryCollision

#endif
