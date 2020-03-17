#include "simulation.hpp"
#include <omp.h>

// Deformation gradient is updated based on the NEW GRID VELOCITIES and the OLD PARTICLE POSITIONS


void Simulation::deformationUpdate_Parallel(){

    unsigned int plastic_count = 0;

    T x0 = grid.x(0);
    T y0 = grid.y(0);

    #pragma omp parallel for reduction(+:plastic_count) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        TM2 sum = TM2::Zero();

        T xp = particles.x(p);
        T yp = particles.y(p);

        unsigned int i_base = std::floor((xp-x0)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((yp-y0)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x(i);
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y(j);

                TV2 vi;
                vi(0) = grid.vx(i,j);
                vi(1) = grid.vy(i,j);

                TV2 grad_wip;
                grad_wip(0) = gradx_wip(xp, yp, xi, yi, one_over_dx);
                grad_wip(1) = grady_wip(xp, yp, xi, yi, one_over_dx);

                sum += vi * grad_wip.transpose();
            } // end loop i
        } // end loop j

        TM2 Fe_trial = particles.F[p];
        Fe_trial = Fe_trial + dt * sum * Fe_trial;
        particles.F[p] = Fe_trial;

        if (plastic_model == VonMises){
            Eigen::JacobiSVD<TM2> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
            TV2 hencky = svd.singularValues().array().log(); // Jixie does not use abs value, however Pradhana-thesis does.
            T   hencky_trace = hencky.sum();
            TV2 hencky_deviatoric = hencky - (hencky_trace / 2.0) * TV2::Ones();
            T   hencky_deviatoric_norm = hencky_deviatoric.norm();

            T delta_gamma = hencky_deviatoric_norm - yield_stress / (2 * mu);
            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.eps_pl_dev(p) += delta_gamma;
                hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        } // end VonMises Plasticity
        else if (plastic_model == NoPlasticity){
            // Do nothing
        }
        else{
            debug("You specified an unvalid PLASTIC model!");
        }
    } // end loop over particles

    debug("               projected particles = ", plastic_count, " / ", Np);

} // end deformationUpdate_Parallel
