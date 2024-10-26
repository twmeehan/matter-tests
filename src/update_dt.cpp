#include "simulation.hpp"

void Simulation::updateDt(){

    // T max_speed = std::sqrt((particles.vx.array().square() + particles.vy.array().square()).maxCoeff());
    auto max_velocity_it = std::max_element( particles.v.begin(), particles.v.end(),
                                             []( const TV &v1, const TV &v2 )
                                             {
                                                 return v1.squaredNorm() < v2.squaredNorm();
                                             } );
    T max_speed = (*max_velocity_it).norm();

    if (max_speed >= wave_speed){
        debug("DETECTED SPEED LARGER THAN ELASTIC WAVE SPEED!!!");
        exit = 1;
        return;
    }

#ifdef WARNINGS
    debug("               dt_max = ", dt_max);
#endif

    if (std::abs(max_speed) > 1e-10){
        T dt_cfl = cfl * dx / max_speed;
#ifdef WARNINGS
        debug("               dt_cfl = ", dt_cfl);
#endif
        dt = std::min(dt_cfl, dt_max);
    } else {
        dt = dt_max;
#ifdef WARNINGS
        debug("               dt_cfl = not computed, max_speed too low");
#endif
    }

    dt = std::min(dt, frame_dt*(frame+1) - time);
    dt = std::min(dt, final_time         - time);

#ifdef WARNINGS
    debug("               dt     = ", dt    );
#endif

    // if (dt > dt_cfl){
    //     debug("TIME STEP IS TOO BIG COMPARED TO CFL!!!");
    //     exit = 1;
    //     return;
    // }
    // if (dt > dt_max){
    //     debug("TIME STEP IS TOO BIG COMPARED TO ELASTIC WAVE SPEED!!!");
    //     exit = 1;
    //     return;
    // }


    if (gravity_special){

        ///////////// ALT 1 /////////////////////////////
        if (time < gravity_time){
            gravity = gravity_final * time/gravity_time;
        }
        else{
            gravity = gravity_final;
            // if (no_liftoff){
            //     for(int p=0; p<Np; p++){
            //         if (particles.x[p](0) > 0.5*Lx){
            //             particles.x[p](0) -= 0.5*Lx;
            //             particles.x[p](1) -= Ly+10*dx;
            //         }
            //     }
            //     no_liftoff = false;
            }
        /////////////////////////////////////////////////

        /////////// ALT 2 ///////////////////////////////
        // T theta_i = 0;
        // T theta_f = 5;
        // T theta = gravity_time * time; // here gravity_time means degrees per sec;
        // theta = std::min(theta, theta_f);
        // debug("theta = ", theta);
        // theta = theta * M_PI / 180;
        // gravity = TV::Zero();
        // gravity[0] = +9.81 * std::sin(theta);
        // gravity[1] = -9.81 * std::cos(theta);
        /////////////////////////////////////////////////

    } // end gravity_special



} // end updateDt
