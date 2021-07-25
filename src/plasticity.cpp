#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DPSimpleSoft || plastic_model == QuadraticLars){

        Eigen::JacobiSVD<TM> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
        // TV hencky = svd.singularValues().array().log(); // VonMises
        TV hencky = svd.singularValues().array().abs().max(1e-4).log(); // QuadraticLars
        T  hencky_trace = hencky.sum();
        TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
        T  hencky_deviatoric_norm = hencky_deviatoric.norm();

        if (hencky_deviatoric_norm > 0)
            hencky_deviatoric /= hencky_deviatoric_norm; // normalize the deviatoric vector so it gives a unit vector specifying the deviatoric direction

        particles.hencky[p] = hencky;

        if (plastic_model == VonMises){

            // Linear Hardening/Softening
            T yield_stress = std::max( (T)1e-3, particles.yield_stress_orig[p] + xi * particles.eps_pl_dev[p] + xi_nonloc * particles.eps_pl_dev_nonloc[p]);

            T delta_gamma = hencky_deviatoric_norm - yield_stress / mu_sqrt6;

            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.delta_gamma[p] = delta_gamma;

                // NB NB NB NB The following 3 lines should be commented out as this is done in plasticity_projection
                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;

            }

        } // end VonMises

        else if (plastic_model == QuadraticLars){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu_sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            T particle_beta;
            T particle_p0_hard;
            if (particles.fail_crit[p]){ // if particle has failed
                particle_beta = 0;
                particle_p0_hard = 100 * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));
            } else{
                if (particles.eps_pl_vol[p] > 0){ //  if particle has not failed in the past, but fails now
                    particle_beta = 0;
                    particle_p0_hard = 100 * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));
                    particles.fail_crit[p] = true;
                } else{ // if particle still elastic
                    particle_beta = beta;
                    particle_p0_hard = p0  * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));
                }
            }

            // bool perform_rma =   CamClayReturnMapping(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, p0_hard, beta, mu, K);
            // bool perform_rma = QuadraticReturnMapping(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, p0_hard, beta, mu, K);
            bool perform_rma = AnalQuadReturnMapping(p_stress, q_stress, exit, M, particle_p0_hard, particle_beta); // p_stress, q_stress will now be the stress at n+1

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T eps_pl_dev_instant = (q_trial - q_stress) / mu_sqrt6;
                particles.eps_pl_dev[p] += eps_pl_dev_instant;

                particles.delta_gamma[p] = eps_pl_dev_instant; // TEMPORARY FOR VIZ.

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p] += eps_pl_vol_inst;

                hencky = q_stress / mu_sqrt6 * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        }

    } // end VonMises or DPSimpleSoft

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

} // end plasticity
