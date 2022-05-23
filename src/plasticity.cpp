#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DruckerPrager || plastic_model == DPSoft || plastic_model == ModifiedCamClay || plastic_model == ModifiedCamClayHard || plastic_model == PerzynaMCC || plastic_model == PerzynaVM || plastic_model == PerzynaDP || plastic_model == PerzynaMuIDP || plastic_model == PerzynaMuIMCC || plastic_model == PerzynaSinterMCC){

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
            // T yield_stress = std::max( (T)1e-3, particles.yield_stress_orig[p] + xi * particles.eps_pl_dev[p] + xi_nonloc * particles.eps_pl_dev_nonloc[p]);

            T yield_stress = std::max( (T)1e-3, yield_stress_orig + xi * particles.eps_pl_dev[p]);

            T delta_gamma = hencky_deviatoric_norm - yield_stress / mu_sqrt6;

            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.delta_gamma[p] = delta_gamma / dt;

                // NB! If using the NONLOCAL approach: The following 3 lines should be commented out as this is done in plasticity_projection
                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;

                // NB sintering
                // particles.sinter_S[p] += dt/sinter_tc*(sinter_Sinf-particles.sinter_S[p]) - particles.sinter_S[p] * delta_gamma / sinter_ec;

            } // end plastic projection

        } // end VonMises

        else if (plastic_model == DruckerPrager){

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

        else if (plastic_model == DPSoft){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = mu_sqrt6 * hencky_deviatoric_norm;

            T p_tip_orig = -dp_cohesion/dp_slope;
            T p_tip      = p_tip_orig * std::exp(-xi * particles.eps_pl_dev[p]);
            T p_shift    = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) <= p_tip){
                plastic_count++;
                T delta_gamma             = hencky_deviatoric_norm;
                particles.delta_gamma[p]  = delta_gamma / dt;
                particles.eps_pl_dev[p]  += delta_gamma;
                T p_proj                  = p_tip_orig * std::exp(-xi * particles.eps_pl_dev[p]); // > p_trial
                T eps_pl_vol_inst         = (p_proj-p_trial)/K;
                particles.eps_pl_vol[p]          += eps_pl_vol_inst;
                particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
            else{ // if right of shifted tip (incl elastic states)
                particles.delta_gamma[p] = 0 / dt; // for the elastic particles, the plastic particles have their delta_gamma overwritten in the next if
            }

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + (-p_tip*dp_slope); // not sure if we should really shift this intersection!!!

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {
                plastic_count++;

                T delta_gamma_0   = hencky_deviatoric_norm - q_yield / mu_sqrt6;
                T temp_eps_pl_dev = particles.eps_pl_dev[p] + delta_gamma_0;
                T p_proj          = p_tip_orig * std::exp(-xi * temp_eps_pl_dev);

                // if left of tip - project from p_trial to p_proj
                if ((p_trial+p_shift) < p_proj){
                    T delta_gamma             = hencky_deviatoric_norm;
                    particles.delta_gamma[p]  = delta_gamma / dt;
                    particles.eps_pl_dev[p]  += delta_gamma;
                    T eps_pl_vol_inst                 = (p_proj-p_trial)/K;
                    particles.eps_pl_vol[p]          += eps_pl_vol_inst;
                    particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
                    hencky = -p_proj/(K*dim) * TV::Ones();
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                }
                // if right of tip - p_trial becomes p
                else{
                    T q_yield_new = dp_slope * (p_trial+p_shift) + (-p_proj*dp_slope) ;
                    T delta_gamma = hencky_deviatoric_norm - q_yield_new / mu_sqrt6;
                    hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                    particles.eps_pl_dev[p] += delta_gamma;
                    particles.delta_gamma[p] = delta_gamma / dt;
                    particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
                }
            } // end plastic projection projection
            // right of tip AND inside yield surface
            else{ // elastic states
                particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
            }

        } // end DPSoft

        else if (plastic_model == PerzynaVM){

            // trial q-stress (in q format)
            T stress = mu_sqrt6 * hencky_deviatoric_norm;

            // update yield stress (q format) based on plastic strain ( first time: yield_stress = yield_stress_orig since epsilon_pl_dev = 0 and exp(..) = 1 )
            T yield_stress = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * particles.eps_pl_dev[p]);

            //// Only for capped von Mises /////
            T p_trial = -K * hencky_trace;
            if (p_trial < vm_ptensile * exp(-xi * particles.eps_pl_vol[p])){
                T delta_gamma = stress / mu_sqrt6;
                T eps_pl_vol_inst = -p_trial/K;
                particles.F[p] = svd.matrixU() * svd.matrixV().transpose();
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma /dt;
            }
            ////////////////////////////////////

            else if (stress > yield_stress) {

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

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;
            } // end plastic projection projection
        } // end PerzynaVM

        else if (plastic_model == PerzynaDP){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = mu_sqrt6 * hencky_deviatoric_norm;

            T p_tip   = -dp_cohesion/dp_slope;
            T p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) < p_tip){
                T delta_gamma     = hencky_deviatoric_norm;
                T p_proj          = p_tip; // > p_trial
                T eps_pl_vol_inst = (p_proj-p_trial)/K;
                plastic_count++;
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.delta_gamma[p]          = delta_gamma;
                particles.eps_pl_dev[p]          += delta_gamma;
                particles.eps_pl_vol[p]          += eps_pl_vol_inst;
                particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
            }
            else{ // if right of shifted tip (incl elastic states)
                particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
                particles.delta_gamma[p] = 0; // for the elastic particles, the plastic particles have their delta_gamma overwritten in the next if
            }

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion; // not sure if we should really shift this intersection!!!

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {

                plastic_count++;

                T delta_gamma = 0.01 * (q_trial - q_yield) / mu_sqrt6; // initial guess

                int max_iter = 60;
                for (int iter = 0; iter < max_iter; iter++) {
                    if (iter == max_iter - 1){ // did not break loop
                        debug("PerzynaDP: FATAL did not exit loop at iter = ", iter);
                        exit = 1;
                    }

                    T tm = perzyna_visc * delta_gamma + dt;
                    T tmp = dt / tm;
                    T tmp1 = std::pow(tmp, perzyna_exp);

                    T residual = (q_trial - mu_sqrt6 * delta_gamma) * tmp1 - q_yield;
                    if (std::abs(residual) < 1e-1) {
                        break;
                    }

                    T residual_diff = -mu_sqrt6 * tmp1 + (q_trial - mu_sqrt6 * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm);

                    if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                        debug("PerzynaDP: residual_diff too small in abs value = ", residual_diff);
                        exit = 1;
                    }

                    delta_gamma -= residual / residual_diff;

                    if (delta_gamma < 0) // not possible and can also lead to division by zero
                        delta_gamma = 1e-10;

                } // end N-R iterations

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;
            } // end plastic projection projection
        } // end PerzynaDP


        else if (plastic_model == PerzynaMuIDP){

            T fac_Q = in_numb_ref * dt / (grain_diameter*std::sqrt(rho_s)); // NB: Use 2 * grain diameter if using the other definiton

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = mu_sqrt6 * hencky_deviatoric_norm;

            T p_tip   = -dp_cohesion/dp_slope;
            T p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) < p_tip){
                T delta_gamma     = hencky_deviatoric_norm;
                T p_proj          = p_tip; // > p_trial
                T eps_pl_vol_inst = (p_proj-p_trial)/K;
                plastic_count++;
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.delta_gamma[p]          = delta_gamma;               // NB!
                particles.eps_pl_dev[p]          += delta_gamma;
                particles.eps_pl_vol[p]          += eps_pl_vol_inst;
                particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
            }
            else{ // if right of shifted tip (incl elastic states)
                particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
                particles.delta_gamma[p] = 0; // for the elastic particles this will be zero, the plastic particles have their delta_gamma overwritten in the next if
            }

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion; // not sure if we should really shift this intersection!!!

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {

                plastic_count++;

                //////////////////////////////////////////////
                // Method 1 - not peric, no mu_i, only exp=1
                //////////////////////////////////////////////

                // T delta_gamma = (q_trial - q_yield) / (perzyna_visc/dt + mu_sqrt6);

                //////////////////////////////////////////////
                // Method 2 - not peric, with mu_i, only exp=1
                //////////////////////////////////////////////
                T fac_a = mu_sqrt6; // always positive
                T fac_b = p_trial*(mu_2-mu_1) + mu_sqrt6*fac_Q*std::sqrt(p_trial) - (q_trial-q_yield);
                T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(p_trial); // always negative

                T delta_gamma = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always psoitive because a>0 and c<0

                // particles.viscosity[p] = 0.5*p_trial*dt*(mu_2-mu_1) / (fac_Q*std::sqrt(p_trial) + delta_gamma);
                particles.viscosity[p] = p_trial*dt*(mu_2-mu_1) / (fac_Q*std::sqrt(p_trial) + delta_gamma);

                // T in_numb = (2*delta_gamma/dt*grain_diameter / std::sqrt(p_trial/rho_s));
                T in_numb = (delta_gamma/dt*grain_diameter / std::sqrt(p_trial/rho_s));
                particles.muI[p] = mu_1 + (mu_2-mu_1) / (in_numb_ref/in_numb + 1);

                //////////////////////////////////////////////
                // Method 2 but linearized mu(I) law
                //////////////////////////////////////////////
                // T visc = (mu_2-mu_1) * grain_diameter * std::sqrt(rho_s*p_trial) / in_numb_ref;
                // T delta_gamma = (q_trial - q_yield) / (visc/dt + mu_sqrt6);
                // particles.viscosity[p] = visc;
                // particles.muI[p] = mu_1 + (mu_2-mu_1) * (2*delta_gamma*grain_diameter / std::sqrt(p_trial/rho_s)) / in_numb_ref;

                //////////////////////////////////////////////
                // Method 3 - not Peric, with mu_i and any exponent. NB: result looks weird for exponent != 1
                //////////////////////////////////////////////

                // // initial guess
                // // T delta_gamma = 0.01 * (q_trial - q_yield) / mu_sqrt6;
                // T delta_gamma = 0.0;
                //
                // int max_iter = 60;
                // for (int iter = 0; iter < max_iter; iter++) {
                //     if (iter == max_iter - 1){ // did not break loop
                //         debug("PerzynaMuIDP: FATAL did not exit loop at iter = ", iter);
                //         // exit = 1;
                //     }
                //     // T tmp = 0.5*p_trial*(mu_2-mu_1) / (fac_Q*std::sqrt(p_trial)/delta_gamma+1);
                //     // T residual = std::pow(tmp, perzyna_exp) - (q_trial - mu_sqrt6*delta_gamma) + q_yield;
                //     // if (std::abs(residual) < 1e-1) {
                //     //     break;
                //     // }
                //     // T d_tmp_d_deltagamma = 0.5 * fac_Q * p_trial * std::sqrt(p_trial) * (mu_2-mu_1) / std::pow(fac_Q * std::sqrt(p) + delta_gamma, 2);
                //     // T residual_diff = perzyna_exp * std::pow(tmp, perzyna_exp-1) * d_tmp_d_deltagamma + mu_sqrt6;
                //
                //     if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                //         debug("PerzynaMuIDP: residual_diff too small in abs value = ", residual_diff);
                //         exit = 1;
                //     }
                //
                //     delta_gamma -= residual / residual_diff;
                //
                //     if (delta_gamma < 0) // not possible and can also lead to division by zero
                //         delta_gamma = 1e-10;
                //
                // } // end N-R iterations

                ////////////////////////////////////////////
                // END METHODS
                ////////////////////////////////////////////

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;

                // // check:
                // T q_new = q_trial - mu_sqrt6 * delta_gamma; // mu_sqrt6 * (hencky - (hencky.sum() / dim) * TV::Ones()).norm();
                // T dg = dt * (q_new-q_yield) / particles.viscosity[p];
                // T rel_error = std::abs(delta_gamma-dg)/dg;
                // if (rel_error > 1e-3 && dg > 1e-13 && delta_gamma > 1e-13){
                //     debug(rel_error, "\t", delta_gamma, "\t", dg);
                //     exit = 1;
                // }

            } // end plastic projection
            else{
                particles.viscosity[p] = 0;
                particles.muI[p] = 0;
            }
        } // end PerzynaMuIDP

        else if (plastic_model == PerzynaMuIMCC) {

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu_sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            // EXLICIT HARDENING
            // T particle_p0 = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]));
            // bool perform_rma = ModifiedCamClayRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K);

            // IMPLICIT HARDENING
            bool perform_rma = ModifiedCamClayHardRMA(p_stress, q_stress, exit, M, particles.eps_pl_vol_mcc[p], beta, mu, K, xi);

            if (perform_rma){
                plastic_count++;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;
                particles.eps_pl_vol_mcc[p] += eps_pl_vol_inst;

                // ONLY if implicit hardening
                T particle_p0 = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]));

                T fac_Q = in_numb_ref / (grain_diameter*std::sqrt(rho_s)); // NB: Use 2 * grain diameter if using the other definiton

                p_stress = std::max(p_stress, -beta*particle_p0);
                p_stress = std::min(p_stress, particle_p0);
                T q_yield = M * std::sqrt( (particle_p0-p_stress)*(beta*particle_p0+p_stress) / (1+2*beta) );

                if (q_trial < q_yield){
                    q_stress = q_yield;
                    debug("Setting q = q_yield = ", q_yield);
                    particles.viscosity[p] = 0;
                    particles.muI[p] = 0;
                }
                else{
                    T fac_a = mu_sqrt6 * dt; // always positive
                    T fac_b = std::abs(p_stress)*(mu_2-mu_1) + mu_sqrt6*dt*fac_Q*std::sqrt(std::abs(p_stress)) - (q_trial-q_yield);
                    T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(std::abs(p_stress)); // always negative

                    T gamma_dot_S = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always positive because a>0 and c<0

                    q_stress = std::max(q_yield, q_trial - mu_sqrt6 * dt * gamma_dot_S);

                    // particles.viscosity[p] = 0.5*std::abs(p_stress)*(mu_2-mu_1) / (fac_Q*std::sqrt(std::abs(p_stress)) + gamma_dot_S);
                    particles.viscosity[p] = std::abs(p_stress)*(mu_2-mu_1) / (fac_Q*std::sqrt(std::abs(p_stress)) + gamma_dot_S);

                    // T in_numb = (2*gamma_dot_S*grain_diameter / std::sqrt(std::abs(p_stress)/rho_s));
                    T in_numb = (gamma_dot_S*grain_diameter / std::sqrt(std::abs(p_stress)/rho_s));
                    particles.muI[p] = mu_1 + (mu_2-mu_1) / (in_numb_ref/in_numb + 1);
                }

                T eps_pl_dev_instant = (q_trial - q_stress) / mu_sqrt6;
                particles.eps_pl_dev[p] += eps_pl_dev_instant;
                particles.delta_gamma[p] = eps_pl_dev_instant / dt;

                hencky = q_stress / mu_sqrt6 * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            } // end plastic projection
            else{
                particles.viscosity[p] = 0;
                particles.muI[p] = 0;
            }


        } // end PerzynaMuIMCC


        else if (plastic_model == ModifiedCamClay || plastic_model == ModifiedCamClayHard || plastic_model == PerzynaMCC || plastic_model == PerzynaSinterMCC){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu_sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            // T particle_p0 = p0;
            // T particle_p0 = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p])); // NB small p0 may be problematic for viscous MCC
            T particle_p0 = std::max(T(1e-2), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]) * (1+particles.sinter_S[p]) );


            ///// HARDNING ALT 3
            // T p0_aftersoft = 20e3;
            // T p0_min = 100;
            // T particle_p0 = p0;
            // T particle_beta = beta;
            // if (particles.eps_pl_vol_abs[p] > 0){ // if plastic
            //     if (particles.fail_crit[p]){ // if finished with softening phase
            //         particle_p0 = std::max(p0_min, (T)(p0_aftersoft * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi))));
            //         particle_beta = 0; // ideally zero
            //     } else { // softening continues
            //         particle_p0 = p0 * std::exp( (1 - std::exp(xi_nonloc*particles.eps_pl_vol_abs[p])) / (rho/1000 * xi) );
            //         particle_beta = beta;
            //         if (particle_p0 < p0_aftersoft){ // softening should stop
            //             particle_p0 = std::max(p0_min, (T)(p0_aftersoft * std::exp((1 - std::exp(particles.eps_pl_vol[p])) / (rho/1000 * xi))));
            //             particle_beta = 0; // ideally zero
            //             particles.fail_crit[p] = true;
            //         }
            //     }
            // } // end if plastic


            bool perform_rma;
            if (plastic_model == ModifiedCamClay)
            {
                perform_rma = ModifiedCamClayRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K);
            }
            else if (plastic_model == ModifiedCamClayHard)
            {
                perform_rma = ModifiedCamClayHardRMA(p_stress, q_stress, exit, M, particles.eps_pl_vol_mcc[p], beta, mu, K, xi);

            }
            else if (plastic_model == PerzynaMCC)
            {
                perform_rma = PerzynaMCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, dt, dim, perzyna_visc);
            }
            else if (plastic_model == PerzynaSinterMCC)
            {
                perform_rma = PerzynaSinterMCCRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, dt, dim, particles.eps_pl_vol[p], particles.sinter_S[p], perzyna_visc, sinter_Sinf, sinter_tc, sinter_ec, xi);
            }

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T eps_pl_dev_instant = (q_trial - q_stress) / mu_sqrt6;
                particles.eps_pl_dev[p] += eps_pl_dev_instant;
                particles.delta_gamma[p] = eps_pl_dev_instant / dt;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;
                particles.eps_pl_vol_mcc[p] += eps_pl_vol_inst;

                T delta_S = dt/sinter_tc*(sinter_Sinf-particles.sinter_S[p]) - particles.sinter_S[p] * std::abs(eps_pl_vol_inst) / sinter_ec;
                particles.sinter_S[p] = std::min(T(sinter_Sinf), std::max(T(0.0), particles.sinter_S[p] + delta_S));

                // particles.eps_pl_vol_abs[p] += std::abs(eps_pl_vol_inst);
                // if (particles.fail_crit[p]){
                //     particles.eps_pl_vol_mcc[p] += -std::abs(eps_pl_vol_inst);
                // }else{
                //     particles.eps_pl_vol_mcc[p] += xi_nonloc * std::abs(eps_pl_vol_inst);
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
