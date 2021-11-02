#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DruckerPrager || plastic_model == Curved || plastic_model == PerzynaVM || plastic_model == PerzynaDP || plastic_model == PerzynaMuIDP){

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
                particles.delta_gamma[p] = delta_gamma;

                // NB! If using the NONLOCAL approach: The following 3 lines should be commented out as this is done in plasticity_projection
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

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;
            } // end plastic projection projection
        } // end PerzynaVM

        else if (plastic_model == PerzynaDP){

            T mu_sqrt6 = mu * 2.44948974278317809819728407471;

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
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion;

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
                        debug("PerzynaDP: residual_diff too small in abs value = ", residual_diff);
                        exit = 1;
                    }

                    delta_gamma -= residual / residual_diff;
                } // end N-R iterations

                hencky -= delta_gamma * hencky_deviatoric; //  note use of delta_gamma instead of delta_gamma_nonloc as in plasticity_projection
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;
            } // end plastic projection projection
        } // end PerzynaDP


        else if (plastic_model == PerzynaMuIDP){

            T mu_sqrt6 = mu * 2.44948974278317809819728407471;

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
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion;

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {

                plastic_count++;

                // T delta_gamma = 0.01 * (q_trial - q_yield) / mu_sqrt6; // initial guess
                T delta_gamma = 0.0;

                int max_iter = 60;
                for (int iter = 0; iter < max_iter; iter++) {
                    if (iter == max_iter - 1){ // did not break loop
                        debug("PerzynaMuIDP: FATAL did not exit loop at iter = ", iter);
                        // exit = 1;
                    }

                    /////////////  Mu(I) Rheology Params  ///////////////
                    // T grain_diameter  = 0.4e-3;
                    // T in_numb_ref     = 2.65;
                    // T mu_1            = 0.38;
                    // T mu_2            = 0.68;
                    // T p_ref           = 1e7;

                    T grain_diameter  = 4;
                    T in_numb_ref     = 2.65;
                    T mu_1            = 0.0;
                    T mu_2            = 0.9;
                    //////////////////////////////////////////////

                    // Method 1

                    // T in_numb = 2 * grain_diameter * delta_gamma / ( dt * std::sqrt(p_trial/rho) );
                    // T mu_i = mu_1 + (mu_2-mu_1) / (in_numb_ref/in_numb+1);
                    // perzyna_visc = mu_i * p_trial * dt / (2 * delta_gamma * p_ref);
                    // // if (perzyna_visc > 10)
                    //     //debug("Visc = ", perzyna_visc);

                    // T tm = perzyna_visc * delta_gamma + dt;
                    // T tmp = dt / tm;
                    // T tmp1 = std::pow(tmp, perzyna_exp);

                    // T residual = (q_trial - mu_sqrt6 * delta_gamma) * tmp1 - q_yield;
                    // if (std::abs(residual) < 1e-1) {
                    //     break;
                    // }
                    //
                    // T d_mui_d_deltagamma  = (mu_2-mu_1) * in_numb_ref * 2*grain_diameter / ( (in_numb_ref+in_numb)*(in_numb_ref+in_numb) * dt * std::sqrt(p_trial/rho) );
                    // T d_visc_d_deltagamma = p_trial * dt * (d_mui_d_deltagamma * delta_gamma - mu_i) / (2 * p_ref * delta_gamma*delta_gamma);
                    // T d_tmp1_d_deltagamma = dt * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * ( perzyna_visc + d_visc_d_deltagamma * delta_gamma ) / ( tm*tm );
                    // T residual_diff = -mu_sqrt6 * tmp1 + (q_trial - mu_sqrt6 * delta_gamma) * d_tmp1_d_deltagamma;

                    // Method 2

                    // T in_numb = 2 * grain_diameter * delta_gamma / ( dt * std::sqrt(p_trial/rho) );
                    // T mu_i = mu_1 + (mu_2-mu_1) / (in_numb_ref/in_numb+1);
                    // perzyna_visc = mu_i * p_trial * dt / (2 * delta_gamma * p_ref); // NOT USED
                    //
                    // T tmp = 2*p_ref / (mu_i*p_trial + 2*p_ref);
                    // T tmp1 = std::pow(tmp, perzyna_exp);
                    //
                    // T residual = (q_trial - mu_sqrt6 * delta_gamma) * tmp1 - q_yield;
                    // if (std::abs(residual) < 1e-1) {
                    //     break;
                    // }
                    //
                    // T half_dimlessp = 0.5*(p_trial/p_ref);
                    // T d_tmp_d_deltagamma = -half_dimlessp * (mu_2-mu_1) / ( (half_dimlessp*mu_i+1) * (half_dimlessp*mu_i+1) ) * in_numb_ref*in_numb / ( (in_numb_ref+in_numb)*(in_numb_ref+in_numb) ) / delta_gamma;
                    // T d_tmp1_d_deltagamma = perzyna_exp * std::pow(tmp, perzyna_exp-1) * d_tmp_d_deltagamma;
                    // T residual_diff = -mu_sqrt6 * tmp1 + (q_trial - mu_sqrt6 * delta_gamma) * d_tmp1_d_deltagamma;

                    // Method 3

                    // T residual = -q_yield + std::pow(dt/(dt + 0.5*dt*p_trial*(mu_1 + (-mu_1 + mu_2)/(1 + 0.5*dt*in_numb_ref*std::sqrt(p_trial)/(delta_gamma*grain_diameter*std::sqrt(rho))))/p_ref), perzyna_exp)*(-delta_gamma*mu_sqrt6 + q_trial);
                    // if (std::abs(residual) < 1e-1) {
                    //     break;
                    // }
                    //
                    // T residual_diff = -mu_sqrt6*std::pow(dt/(dt + 0.5*dt*p_trial*(mu_1 + (-mu_1 + mu_2)/(1 + 0.5*dt*in_numb_ref*std::sqrt(p_trial)/(delta_gamma*grain_diameter*std::sqrt(rho))))/p_ref), perzyna_exp) - 1.0L/4.0L*std::pow(dt, 2)*in_numb_ref*std::pow(p_trial, 3.0L/2.0L)*perzyna_exp*std::pow(dt/(dt + 0.5*dt*p_trial*(mu_1 + (-mu_1 + mu_2)/(1 + 0.5*dt*in_numb_ref*std::sqrt(p_trial)/(delta_gamma*grain_diameter*std::sqrt(rho))))/p_ref), perzyna_exp)*(-mu_1 + mu_2)*(-delta_gamma*mu_sqrt6 + q_trial)/(std::pow(delta_gamma, 2)*grain_diameter*p_ref*std::sqrt(rho)*std::pow(1 + 0.5*dt*in_numb_ref*std::sqrt(p_trial)/(delta_gamma*grain_diameter*std::sqrt(rho)), 2)*(dt + 0.5*dt*p_trial*(mu_1 + (-mu_1 + mu_2)/(1 + 0.5*dt*in_numb_ref*std::sqrt(p_trial)/(delta_gamma*grain_diameter*std::sqrt(rho))))/p_ref));


                    // Method 4 - not peric, no mu_i, only exp=1

                    // delta_gamma = (q_trial - q_yield) / (perzyna_visc/dt + mu_sqrt6);
                    // T residual = 0;
                    // T residual_diff = 1;
                    // break;

                    // Method 5 - not peric, with mu_i, only exp=1

                    T fac_Q = in_numb_ref * dt / (2*grain_diameter*std::sqrt(rho));

                    T fac_a = mu_sqrt6;
                    // T fac_b = 0.5*p_trial*(mu_2-mu_1) + mu_sqrt6*fac_Q*std::sqrt(p_trial) - (q_trial-q_yield-0.5*p_trial*mu_1);
                    T fac_b = 0.5*mu_2*p_trial + mu_sqrt6*fac_Q*std::sqrt(p_trial) - (q_trial-q_yield);
                    T fac_c = -(q_trial-q_yield-0.5*p_trial*mu_1) * fac_Q * std::sqrt(p_trial);

                    T delta_gamma_neg = (-fac_b - std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a);
                    T delta_gamma_pos = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a);

                    delta_gamma = std::max(delta_gamma_pos, delta_gamma_neg);
                    if (delta_gamma < 0){
                        // debug("PerzynaMuIDP: delta_gamma = ", delta_gamma);
                        // exit = 1;
                        delta_gamma = 0;
                    }

                    T residual = 0;
                    T residual_diff = 1;
                    break;

                    // Method 6 - not Peric

                    // T in_numb = 2 * grain_diameter * delta_gamma / ( dt * std::sqrt(p_trial/rho) );
                    // T mu_i = mu_1 + in_numb * (mu_2-mu_1) / (in_numb_ref+in_numb);
                    //
                    // T tmp  = 0.5*p_trial*mu_i; // = visc * delta_gamma / dt
                    // T residual = std::pow(tmp, perzyna_exp) - q_trial + mu_sqrt6*delta_gamma + q_yield;
                    //
                    // if (std::abs(residual) < 1e-1) {
                    //     break;
                    // }
                    //
                    // T upper = p_trial*(mu_2-mu_1) * in_numb_ref * dt * std::sqrt(p_trial/rho);
                    // T lower = 4*grain_diameter * std::pow( in_numb_ref * dt * std::sqrt(p_trial/rho) / (2*grain_diameter) + delta_gamma, 2);
                    // T d_tmp_d_deltagamma = upper / lower;
                    // T residual_diff = perzyna_exp * std::pow(tmp, perzyna_exp-1) * d_tmp_d_deltagamma + mu_sqrt6;
                    //
                    // perzyna_visc = tmp*dt/delta_gamma; // for viz. only


                    ////  End Methods ////

                    if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                        debug("PerzynaMuIDP: residual_diff too small in abs value = ", residual_diff);
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
        } // end PerzynaMuIDP


        else if (plastic_model == Curved){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu_sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            T particle_beta = beta;
            // T particle_p0_hard = p0;
            T particle_p0_hard = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p])); // NB small p0 may be problematic for viscous MCC

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


            // bool perform_rma = PerzynaCamClayRMA(p_stress, q_stress, exit, M, particle_p0_hard, particle_beta, mu, K, dt, dim, perzyna_visc);
            // bool perform_rma = PerzynaQuadRMA(p_stress, q_stress, exit, M, particle_p0_hard, particle_beta, mu, K, dt, dim, perzyna_visc);
            bool perform_rma = CamClayRMA(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, particle_p0_hard, particle_beta, mu, K);
            // bool perform_rma = QuadRMA(p_stress, q_stress, exit, hencky_trace, hencky_deviatoric_norm, M, p0_hard, beta, mu, K);
            // bool perform_rma = QuadAnalyticRMA(p_stress, q_stress, exit, M, particle_p0_hard, particle_beta); // p_stress, q_stress will now be the stress at n+1

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T eps_pl_dev_instant = (q_trial - q_stress) / mu_sqrt6;
                particles.eps_pl_dev[p] += eps_pl_dev_instant;
                particles.delta_gamma[p] = eps_pl_dev_instant; // TEMPORARY FOR VIZ.

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;
                particles.eps_pl_vol_mcc[p] += eps_pl_vol_inst;
                // particles.eps_pl_vol_abs[p] += std::abs(eps_pl_vol_inst);
                //
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
