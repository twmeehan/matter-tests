#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DruckerPrager || plastic_model == Curved || plastic_model == PerzynaVM || plastic_model == PerzynaNA){

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

            // NB: q-format. Linear Hardening/Softening
            T yield_stress = std::max( (T)1e-3, particles.yield_stress_orig[p] + xi * particles.eps_pl_dev[p] + xi_nonloc * particles.eps_pl_dev_nonloc[p]);

            T delta_gamma = hencky_deviatoric_norm - yield_stress / mu_sqrt6;

            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.delta_gamma[p] = delta_gamma;

                // NB NB NB NB The following 3 lines should be commented out as this is done in plasticity_projection
                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;

            } // end plastic projection

        } // end VonMises

        else if (plastic_model == DruckerPrager){

            T mu_sqrt6 = mu * 2.44948974278317809819728407471;

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = mu_sqrt6 * hencky_deviatoric_norm;

            T q_yield = dp_slope * p_trial + dp_cohesion;

            // if left of tip
            if (q_yield < 1e-10){
                T delta_gamma = hencky_deviatoric_norm;
                T p_proj = -dp_cohesion/dp_slope; // larger than p_trial
                plastic_count++;
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.delta_gamma[p] = delta_gamma;
                particles.eps_pl_dev[p] += delta_gamma;
                particles.eps_pl_vol[p] += (p_proj-p_trial)/K;
            }
            else{ // right of tipe
                T delta_gamma = hencky_deviatoric_norm - q_yield / mu_sqrt6;

                if (delta_gamma > 0){ // project to yield surface
                    plastic_count++;
                    particles.delta_gamma[p] = delta_gamma;

                    // NB NB NB NB The following 3 lines should be commented out as this is done in plasticity_projection
                    hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                    particles.eps_pl_dev[p] += delta_gamma;
                }
            } // if else side of tip

        } // end DruckerPrager

        else if (plastic_model == PerzynaVM){

            T mu_sqrt6 = mu * 2.44948974278317809819728407471;

            // trial q-stress (in q format)
            T stress = mu_sqrt6 * hencky_deviatoric_norm;

            // update yield stress (q format) based on plastic strain ( first time: yield_stress = yield_stress_orig since epsilon_pl_dev = 0 and exp(..) = 1 )
            T yield_stress = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * particles.eps_pl_dev[p]);

            if (stress > yield_stress) {

                plastic_count++;

                T delta_gamma = 0.1 * (stress - yield_stress) / mu_sqrt6; // initial guess

                int max_iter = 60;
                for (int iter = 0; iter < max_iter; iter++) {
                    if (iter == max_iter - 1){ // did not break loop
                        debug("PerzynaVM: FATAL did not exit loop at iter = ", iter);
                        exit = 1;
                    }

                    if (delta_gamma < 0) // not possible and can also lead to division by zero
                        delta_gamma = 1e-10;

                    T tm = perzyna_visc * delta_gamma + dt;
                    T tmp = dt / tm;
                    T tmp1 = std::pow(tmp, perzyna_exp);

                    T yield_stress_new = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * (particles.eps_pl_dev[p] + delta_gamma));

                    T residual = (stress - mu_sqrt6 * delta_gamma) * tmp1 - yield_stress_new;
                    if (std::abs(residual) < 1e-1) {
                        break;
                    }

                    T yield_stress_new_diff = -xi * (yield_stress_orig - yield_stress_min) * exp(-xi * (particles.eps_pl_dev[p] + delta_gamma));
                    T residual_diff         = -mu_sqrt6 * tmp1 + (stress - mu_sqrt6 * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm) - yield_stress_new_diff;

                    if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                        debug("PerzynaVM: residual_diff too small in abs value = ", residual_diff);
                        exit = 1;
                    }

                    delta_gamma -= residual / residual_diff;
                } // end N-R iterations

                particles.delta_gamma[p] = delta_gamma;

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
            } // end plastic projection projection
        } // end PerzynaVM

        else if (plastic_model == PerzynaNA){

            T mu_sqrt6 = mu * 2.44948974278317809819728407471;

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = mu_sqrt6 * hencky_deviatoric_norm;

            // If DP
            T q_yield = std::max((T)1e-10, dp_slope * p_trial + dp_cohesion);

            // if Non Ass MCC
            // T q_yield;
            // if (p_trial <= -beta*p0){
            //     q_yield = 1e-10;
            // } else if (p_trial >= p0){
            //     q_yield = 1e-10;
            // } else{
            //     q_yield = M*std::sqrt( (p0-p_trial)*(beta*p0+p_trial) / (1+2*beta) );
            // }

            if (q_trial > q_yield) {

                plastic_count++;

                T delta_gamma = 0.1 * (q_trial - q_yield) / mu_sqrt6; // initial guess

                int max_iter = 60;
                for (int iter = 0; iter < max_iter; iter++) {
                    if (iter == max_iter - 1){ // did not break loop
                        debug("PerzynaNA: FATAL did not exit loop at iter = ", iter);
                        exit = 1;
                    }

                    if (delta_gamma < 0) // not possible and can also lead to division by zero
                        delta_gamma = 1e-10;

                    T tm = perzyna_visc * delta_gamma + dt;
                    T tmp = dt / tm;
                    T tmp1 = std::pow(tmp, perzyna_exp);

                    T residual = (q_trial - mu_sqrt6 * delta_gamma) * tmp1 - q_yield;
                    if (std::abs(residual) < 1e-1) {
                        break;
                    }

                    T residual_diff = -mu_sqrt6 * tmp1 + (q_trial - mu_sqrt6 * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm);

                    if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                        debug("PerzynaNA: residual_diff too small in abs value = ", residual_diff);
                        exit = 1;
                    }

                    delta_gamma -= residual / residual_diff;
                } // end N-R iterations

                particles.delta_gamma[p] = delta_gamma;

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
            } // end plastic projection projection
        } // end PerzynaDP

        else if (plastic_model == Curved){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu_sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            T particle_beta = beta;
            T particle_p0_hard = p0;

            ////// HARDNING ALT 0
            // T particle_beta = beta;
            // T particle_p0_hard = p0 * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));

            ////// HARDNING ALT 1
            // T particle_beta = beta;
            // T particle_p0_hard = p0;
            // if (particles.fail_crit[p]){ // if particle has failed
            //     particle_beta = 0;
            //     particle_p0_hard = p0/2 * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));
            // } else{
            //     if (particles.eps_pl_dev[p] > 0){ //  if particle has not failed in the past, but fails now
            //         particle_beta = 0;
            //         particle_p0_hard = p0/2 * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi));
            //         particles.fail_crit[p] = true;
            //     }
            // }

            ////// HARDNING ALT 2
            // T particle_p0_hard = std::max(T(100), (T)(p0  * std::exp(0.0001 * xi * particles.eps_pl_vol[p])));
            // T particle_pt_hard =              (beta * p0) * std::exp(      xi * particles.eps_pl_vol_abs[p]);
            // T particle_beta = particle_pt_hard / particle_p0_hard;

            ///// HARDNING ALT 3
            // T p0_aftersoft = 1000;
            // T p0_min = 1000;
            // T particle_p0_hard = p0;
            // T particle_beta = beta;
            // if (particles.eps_pl_vol_2[p] > 0){ // if plastic
            //     if (particles.fail_crit[p]){ // if finished with softening phase
            //         particle_p0_hard = std::max(p0_min, (T)(p0_aftersoft * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi))));
            //         particle_beta = 0;
            //     } else { // softening continues
            //         particle_p0_hard = p0 * std::exp( (1 - std::exp(xi_nonloc*particles.eps_pl_vol_2[p])) / (rho/1000 * xi) );
            //         particle_beta = beta;
            //         if (particle_p0_hard < p0_aftersoft){ // softening should stop
            //             particle_p0_hard = std::max(p0_min, (T)(p0_aftersoft * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi))));
            //             particle_beta = 0;
            //             particles.fail_crit[p] = true;
            //         }
            //     }
            // } // end if plastic

            ///// HARDNING ALT MCC
            // T p0_aftersoft = 1000.0;
            // T eps_pl_vol_limit = -std::asinh(p0_aftersoft/K) / xi;
            // T particle_p0_hard = K*std::sinh(xi*std::max(-particles.eps_pl_vol_3[p], -eps_pl_vol_limit));
            // T particle_beta = beta;
            // if (particles.eps_pl_dev[p] > 0){ // if plastic
            //     if (particles.fail_crit[p]){ // if finished with softening phase
            //         particle_beta = 0;
            //     } else { // softening continues
            //         if (particle_p0_hard < (p0_aftersoft+1e-3)){ // softening should stop
            //             particle_beta = 0;
            //             particle_p0_hard = p0_aftersoft;
            //             particles.eps_pl_vol_3[p] = eps_pl_vol_limit;
            //             particles.fail_crit[p] = true;
            //         }
            //     }
            // } // end if plastic

            ///// HARDNING ALT MCC 2
            // T p0_aftersoft = 1000;
            //
            // // if elastic or softening phase
            // T particle_p0_hard = K*std::sinh(-xi*particles.eps_pl_vol_3[p]); // initially eps_pl_vol_3 must be negative!!!
            // T particle_beta = beta;
            //
            // // if plastic
            // if (particles.eps_pl_dev[p] > 0){
            //     if (particles.fail_crit[p]){ // if finished with softening phase
            //         particle_beta = 0;
            //         particle_p0_hard = std::max(p0_aftersoft, K*std::sinh(-xi*particles.eps_pl_vol_3[p]));
            //     } else { // softening continues
            //         if (particle_p0_hard < p0_aftersoft){ // softening should stop
            //             particle_beta = 0;
            //             particle_p0_hard = p0_aftersoft;
            //             particles.eps_pl_vol_3[p] = -std::asinh(p0_aftersoft/K) / xi;;
            //             particles.fail_crit[p] = true;
            //         }
            //     }
            // } // end if plastic


            bool perform_rma =   PerzynaQuadReturnMapping(p_stress, q_stress, exit, M, p0, beta, mu, K, dt, dim, perzyna_visc);
            // bool perform_rma = CamClayReturnMapping(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, particle_p0_hard, particle_beta, mu, K);
            // bool perform_rma = QuadraticReturnMapping(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, p0_hard, beta, mu, K);
            // bool perform_rma = AnalQuadReturnMapping(p_stress, q_stress, exit, M, particle_p0_hard, particle_beta); // p_stress, q_stress will now be the stress at n+1

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T eps_pl_dev_instant = (q_trial - q_stress) / mu_sqrt6;
                particles.eps_pl_dev[p] += eps_pl_dev_instant;

                particles.delta_gamma[p] = eps_pl_dev_instant; // TEMPORARY FOR VIZ.

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                // particles.eps_pl_vol_2[p] += std::abs(eps_pl_vol_inst);
                //
                // if (particles.fail_crit[p]){
                //     particles.eps_pl_vol_3[p] += -std::abs(eps_pl_vol_inst);
                // }else{
                //     particles.eps_pl_vol_3[p] += xi_nonloc * std::abs(eps_pl_vol_inst);
                // }

                hencky = q_stress / mu_sqrt6 * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        }

    } // end plastic_model type

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

} // end plasticity
