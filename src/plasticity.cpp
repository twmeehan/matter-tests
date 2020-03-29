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

            // Softening
            T yield_stress = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * particles.eps_pl_dev[p]);

            /////////// REGULARIZATION ///////////
            // T yield_stress_reg =                yield_stress - reg_const_length_sq * particles.regularization(p);
            T yield_stress_reg = std::max( (T)0.0, yield_stress - reg_const_length_sq * particles.regularization[p] );
            if (yield_stress_reg < 0.0){
                debug("NEGATIVE YIELD STRESS!!!");
                exit = 1;
                return;
            }
            //////////////////////////////////////

            T delta_gamma = hencky_deviatoric_norm - yield_stress_reg / (2 * mu);
            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
            }
        } // end VonMises

        else if (plastic_model == DPSimpleSoft){

            T cohesion_proj = particles.cohesion_proj[p];
            if (hencky_trace >= dim * cohesion_proj) { // Project to tip
                plastic_count++;
                cohesion_proj = cohesion * std::exp(-xi * particles.eps_pl_dev[p]);
                particles.F[p] = svd.matrixU() * std::exp(cohesion_proj) * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += hencky_deviatoric_norm;
                particles.eps_pl_vol[p] += (hencky_trace - dim * cohesion_proj);
            }
            else{ // Right of tip
                T delta_gamma = hencky_deviatoric_norm + alpha_K_d_over_2mu * hencky_trace - alpha_K_d_over_2mu * dim * cohesion_proj;
                if (delta_gamma > 0) { // outside yield surface
                    plastic_count++;
                    cohesion_proj = cohesion * std::exp(-xi * particles.eps_pl_dev[p]);
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
