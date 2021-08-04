#include "simulation.hpp"

#ifdef THREEDIM

void Simulation::boundaryCollision(T xi, T yi, T zi, TV& vi){

    // TODO: Make this a vector-based function

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
                    // normal component (y) must be set to zero
                    vy_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::sqrt(vx_rel*vx_rel + vz_rel*vz_rel) < friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vy_rel);
                        else
                            vx_rel += friction * std::abs(vy_rel);
                        if (vz_rel > 0)
                            vz_rel -= friction * std::abs(vy_rel);
                        else
                            vz_rel += friction * std::abs(vy_rel);
                    }
                } // end top and bottom plate
                else if (obj.plate_type == left || obj.plate_type == right){
                    // tangential velocity is the (y,z) components
                    // normal component (x) must be set to zero
                    vx_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::sqrt(vy_rel*vy_rel + vz_rel*vz_rel) < friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vx_rel);
                        else
                            vy_rel += friction * std::abs(vx_rel);
                        if (vz_rel > 0)
                            vz_rel -= friction * std::abs(vx_rel);
                        else
                            vz_rel += friction * std::abs(vx_rel);
                    }
                } // end left or right plate
                else if (obj.plate_type == front || obj.plate_type == back){
                    // tangential velocity is the (x,y) components
                    // normal component (z) must be set to zero
                    vz_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel) < friction * std::abs(vx_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                        vy_rel = 0;
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vz_rel);
                        else
                            vx_rel += friction * std::abs(vz_rel);
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vz_rel);
                        else
                            vy_rel += friction * std::abs(vz_rel);
                    }
                } // end front or back plate
            } // end SLIP

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
                    // tangential velocity is the (x) component
                    // normal component (y) must be set to zero
                    vy_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vy_rel);
                        else
                            vx_rel += friction * std::abs(vy_rel);
                    }
                } // end top and bottom plate
                else if (obj.plate_type == left || obj.plate_type == right){
                    // tangential velocity is the (y) component
                    // normal component (x) must be set to zero
                    vx_rel = 0;
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vx_rel);
                        else
                            vy_rel += friction * std::abs(vx_rel);
                    }
                } // end left or right plate
            } // end SLIP

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
