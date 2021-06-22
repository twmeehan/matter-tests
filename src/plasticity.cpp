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

        if (plastic_model == VonMises){

            // Exponential Softening
            // T yield_stress_local = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * particles.eps_pl_dev[p]);

            // Linear Softening
            T yield_stress_local = yield_stress_orig + xi * particles.eps_pl_dev[p];

            /////////// REGULARIZATION ///////////
            T yield_stress = std::max( (T)1e-3, yield_stress_local + l_sq * particles.reg_laplacian[p] );
            // if (yield_stress_reg < 0.0){
            //     debug("NEGATIVE YIELD STRESS!!!");
            //     exit = 1;
            //     return;
            // }
            //////////////////////////////////////

            T delta_gamma = hencky_deviatoric_norm - yield_stress / (2 * mu);
            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev_inst[p]  = delta_gamma;
                particles.eps_pl_dev[p]      += delta_gamma;
            }

            particles.reg_variable[p] = 1.0;
        } // end VonMises

        else if (plastic_model == DPSimpleSoft){

            if (hencky_trace >= dim * particles.cohesion_proj[p]) { // Project to tip
                plastic_count++;
                particles.cohesion_proj[p] = cohesion * std::exp(-xi * particles.eps_pl_dev[p]);
                particles.F[p] = svd.matrixU() * std::exp(particles.cohesion_proj[p]) * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += hencky_deviatoric_norm;
                particles.eps_pl_vol[p] += (hencky_trace - dim * particles.cohesion_proj[p]);
            }
            else{ // Right of tip
                T delta_gamma = hencky_deviatoric_norm + alpha_K_d_over_2mu * hencky_trace - alpha_K_d_over_2mu * dim * particles.cohesion_proj[p];
                if (delta_gamma > 0) { // outside yield surface
                    plastic_count++;
                    particles.cohesion_proj[p] = cohesion * std::exp(-xi * particles.eps_pl_dev[p]);
                    hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                    particles.eps_pl_dev[p] += delta_gamma;
                } // end outside yield surface
            } // end left or right of tip

        } // end DPSimpleSoft

    } // end VonMises or DPSimpleSoft

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

    // particles.eps_pl_dev[p] = 0.1; // Testing regularization`

} // end plasticity
