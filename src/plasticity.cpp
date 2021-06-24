#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DPSimpleSoft){
        Eigen::JacobiSVD<TM> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
        TV hencky = svd.singularValues().array().log(); // Jixie does not use abs value, however Pradhana-thesis does.
        T  hencky_trace = hencky.sum();
        TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
        T  hencky_deviatoric_norm = hencky_deviatoric.norm();

        particles.hencky[p] = hencky;

        if (plastic_model == VonMises){

            // Linear Hardening/Softening
            T yield_stress = yield_stress_orig + xi * particles.eps_pl_dev[p];

            T delta_gamma = hencky_deviatoric_norm - yield_stress / (2 * mu);
            particles.delta_gamma[p] = delta_gamma;

            if (delta_gamma > 0) // project to yield surface
                plastic_count++;

        } // end VonMises

    } // end VonMises or DPSimpleSoft

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

} // end plasticity
