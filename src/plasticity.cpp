#include "simulation.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VonMises || plastic_model == DruckerPrager || plastic_model == DPSoft || plastic_model == MCC || plastic_model == MCCHard || plastic_model == MCCHardExp || plastic_model == PerzynaMCC || plastic_model == PerzynaVM || plastic_model == PerzynaDP || plastic_model == PerzynaMuIDP || plastic_model == PerzynaMuIMCC  || plastic_model == PerzynaMCCHard || plastic_model == SinterMCC){

        Eigen::JacobiSVD<TM> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
        // TV hencky = svd.singularValues().array().log(); // VonMises
        TV hencky = svd.singularValues().array().abs().max(1e-4).log(); // QuadraticLars
        T  hencky_trace = hencky.sum();
        TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
        T  hencky_deviatoric_norm = hencky_deviatoric.norm();

        if (hencky_deviatoric_norm > 0)
            hencky_deviatoric /= hencky_deviatoric_norm; // normalize the deviatoric vector so it gives a unit vector specifying the deviatoric direction

        // particles.hencky[p] = hencky;

        if (plastic_model == VonMises){

            // NB: q-format. Linear Hardening/Softening
            // T yield_stress = std::max( (T)1e-3, particles.yield_stress_orig[p] + xi * particles.eps_pl_dev[p] + xi_nonloc * particles.eps_pl_dev_nonloc[p]);

            T yield_stress = std::max( (T)1e-3, yield_stress_orig + xi * particles.eps_pl_dev[p]);

            T delta_gamma = hencky_deviatoric_norm - yield_stress / (mu*sqrt6);

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
            T q_trial = mu*sqrt6 * hencky_deviatoric_norm;

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
                T delta_gamma = hencky_deviatoric_norm - q_yield / (mu*sqrt6);

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
            T q_trial = mu*sqrt6 * hencky_deviatoric_norm;

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

                T delta_gamma_0   = hencky_deviatoric_norm - q_yield / (mu*sqrt6);
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
                    T delta_gamma = hencky_deviatoric_norm - q_yield_new / (mu*sqrt6);
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
            T stress = mu*sqrt6 * hencky_deviatoric_norm;

            // update yield stress (q format) based on plastic strain ( first time: yield_stress = yield_stress_orig since epsilon_pl_dev = 0 and exp(..) = 1 )
            T yield_stress = yield_stress_min + (yield_stress_orig - yield_stress_min) * exp(-xi * particles.eps_pl_dev[p]);

            //// Only for capped von Mises /////
            T p_trial = -K * hencky_trace;
            if (p_trial < vm_ptensile * exp(-xi * particles.eps_pl_vol[p])){
                T delta_gamma = stress / (mu*sqrt6);
                T eps_pl_vol_inst = -p_trial/K;
                particles.F[p] = svd.matrixU() * svd.matrixV().transpose();
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma /dt;
            }
            ////////////////////////////////////

            else if (stress > yield_stress) {

                plastic_count++;

                T delta_gamma = 0.1 * (stress - yield_stress) / (mu*sqrt6); // initial guess

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

                    T residual = (stress - mu*sqrt6 * delta_gamma) * tmp1 - yield_stress_new;
                    if (std::abs(residual) < 1e-1) {
                        break;
                    }

                    T yield_stress_new_diff = -xi * (yield_stress_orig - yield_stress_min) * exp(-xi * (particles.eps_pl_dev[p] + delta_gamma));
                    T residual_diff         = -mu*sqrt6 * tmp1 + (stress - mu*sqrt6 * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm) - yield_stress_new_diff;

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
            T q_trial = mu*sqrt6 * hencky_deviatoric_norm;

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

                T delta_gamma = 0.01 * (q_trial - q_yield) / (mu*sqrt6); // initial guess

                int max_iter = 60;
                for (int iter = 0; iter < max_iter; iter++) {
                    if (iter == max_iter - 1){ // did not break loop
                        debug("PerzynaDP: FATAL did not exit loop at iter = ", iter);
                        exit = 1;
                    }

                    T tm = perzyna_visc * delta_gamma + dt;
                    T tmp = dt / tm;
                    T tmp1 = std::pow(tmp, perzyna_exp);

                    T residual = (q_trial - mu*sqrt6 * delta_gamma) * tmp1 - q_yield;
                    if (std::abs(residual) < 1e-1) {
                        break;
                    }

                    T residual_diff = -mu*sqrt6 * tmp1 + (q_trial - mu*sqrt6 * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm);

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


        // else if ((plastic_model == PerzynaMuIMCC) && (particles.x[p](0) < 1.2*Lx)){
        else if (plastic_model == PerzynaMuIDP){

            T e_mu_prefac, g_mu_prefac, dg_prefac;
            if (use_jop_definitions){
                e_mu_prefac = sqrt2 * mu;
                // g_mu_prefac = 2.0 * mu;
                // dg_prefac = 1.0 / sqrt2;
                g_mu_prefac =  mu;
                dg_prefac = sqrt2;
            } else{
                e_mu_prefac = sqrt6 * mu;
                g_mu_prefac = e_mu_prefac;
                dg_prefac = 1;
            }

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = e_mu_prefac * hencky_deviatoric_norm;

            T p_tip   = -dp_cohesion/dp_slope;
            T p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            particles.muI[p]         = mu_1;
            particles.viscosity[p]   = 0;
            particles.delta_gamma[p] = 0;

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) < p_tip){
                T delta_gamma     = dg_prefac * hencky_deviatoric_norm;
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
            }

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion; // not sure if we should really shift this intersection!!!

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {

                plastic_count++;

                T fac_a = g_mu_prefac * dt; // always positive
                T fac_b = p_trial*(mu_2-mu_1) + g_mu_prefac*dt*fac_Q*std::sqrt(std::abs(p_trial)) - (q_trial-q_yield);
                T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(std::abs(p_trial)); // always negative

                // this is gamma_dot:
                T delta_gamma = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always positive because a>0 and c<0

                T mu_i                 = mu_1 + (mu_2 - mu_1) / (fac_Q * std::sqrt(std::abs(p_trial)) / delta_gamma + 1.0);
                particles.muI[p]       = mu_i;
                particles.viscosity[p] = (mu_i - mu_1) * std::abs(p_trial) / delta_gamma;

                delta_gamma *= dt; // this is the actual delta_gamma

                hencky -= (1.0/dg_prefac) * delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;
                particles.delta_gamma[p] = delta_gamma;

            } // end plastic projection

        } // end PerzynaMuIDP

        else if (plastic_model == PerzynaMuIMCC) {

            T e_mu_prefac, g_mu_prefac, dg_prefac, rma_prefac;
            if (use_jop_definitions){
                e_mu_prefac = sqrt2 * mu;
                // g_mu_prefac =   2.0 * mu; // Jop
                g_mu_prefac = mu; // Kamrin
                rma_prefac = 1;
            } else{
                e_mu_prefac = sqrt6 * mu;
                g_mu_prefac = e_mu_prefac;
                rma_prefac = 3;
            }

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = e_mu_prefac * hencky_deviatoric_norm;

            // make copies
            T q_trial = q_stress;

            // T particle_xi = xi;
            // if (particles.x[p](0) < 1.001*Lx)
            //     particle_xi = xi*10;

            //////////////////////////////////////////////////////////////////////
            /////// IMPLICIT HARDENING
            // bool perform_rma = MCCHardRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, rma_prefac, particles.eps_pl_vol[p]);
            bool perform_rma = MCCHardExpRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, rma_prefac, particles.eps_pl_vol[p]);
            /////// EXLICIT HARDENING
            // T particle_p0 = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]));
            // T particle_p0 = std::max(T(1e-3), K*std::sinh(-xi*particles.eps_pl_vol[p] + std::asinh(p0/K)));
            // T particle_p0 = std::max(T(1e-3), p0*std::sinh(-xi*particles.eps_pl_vol[p] + std::asinh(1.0)));
            // T particle_p0 = std::max(T(1e-2), p0*std::exp(-xi*particles.eps_pl_vol[p]));
            // T particle_p0 = std::max( T(1e-2), p0*std::exp(xi*(1-std::exp(particles.eps_pl_vol[p]))) );
            // T particle_p0 = std::max(T(1e-2), (particles.eps_pl_vol[p] < 0) ? p0*(1.0-std::sinh(xi*particles.eps_pl_vol[p])) : p0*(1.0-std::tanh(xi*particles.eps_pl_vol[p])) );
            // T particle_p0 = std::max(T(1e-3), (particles.eps_pl_vol_mcc[p] < 0) ? K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]) : K*std::tanh(-xi*particles.eps_pl_vol_mcc[p]) );
            // bool perform_rma = MCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, rma_prefac);
            //////////////////////////////////////////////////////////////////////

            particles.muI[p] = mu_1;
            particles.viscosity[p] = 0;

            if (perform_rma){
                plastic_count++;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;
                particles.eps_pl_vol_mcc[p] += eps_pl_vol_inst;

                T q_yield = q_stress;
                if (q_trial < (q_yield + 1e-5)) {
                    q_stress = q_yield;
                    // debug("Setting q = q_yield = ", q_yield);
                }
                else{
                    T fac_a = g_mu_prefac * dt; // always positive
                    T fac_b = std::abs(p_stress)*(mu_2-mu_1) + g_mu_prefac*dt*fac_Q*std::sqrt(std::abs(p_stress)) - (q_trial-q_yield);
                    T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(std::abs(p_stress)); // always negative

                    T gamma_dot_S = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always positive because a>0 and c<0

                    q_stress = std::max(q_yield, q_trial - g_mu_prefac * dt * gamma_dot_S);

                    T mu_i                 = mu_1 + (mu_2 - mu_1) / (fac_Q * std::sqrt(std::abs(p_stress)) / gamma_dot_S + 1.0);
                    particles.muI[p]       = mu_i;
                    particles.viscosity[p] = (mu_i - mu_1) * std::abs(p_stress) / gamma_dot_S;
                }

                T dg_instant = (q_trial - q_stress) / g_mu_prefac;
                particles.eps_pl_dev[p] += dg_instant;
                particles.delta_gamma[p] = dg_instant / dt;

                hencky = q_stress / e_mu_prefac * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            } // end plastic projection

        } // end PerzynaMuIMCC


        else if (plastic_model == MCC || plastic_model == MCCHard || plastic_model == MCCHardExp || plastic_model == PerzynaMCC || plastic_model == SinterMCC || plastic_model == PerzynaMCCHard){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = mu*sqrt6 * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            // T particle_p0 = p0;
            // T particle_p0 = std::max(T(1e-2), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]));
            T particle_p0 = std::max(T(1e-2), p0*std::exp(-xi*particles.eps_pl_vol[p]));
            // T particle_p0 = std::max(T(1e-2), K*std::sinh(-xi*particles.eps_pl_vol_mcc[p]) * (1+particles.sinter_S[p]) );

            bool perform_rma;
            if (plastic_model == MCC)
            {
                perform_rma = MCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, 1);
            }
            else if (plastic_model == MCCHard)
            {
                perform_rma = MCCHardRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, 1, particles.eps_pl_vol[p]);
            }
            else if (plastic_model == MCCHardExp)
            {
                perform_rma = MCCHardExpRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, 1, particles.eps_pl_vol[p]);
            }
            else if (plastic_model == PerzynaMCC)
            {
                perform_rma = PerzynaMCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, dt, dim, perzyna_visc);
            }
            else if (plastic_model == PerzynaMCCHard)
            {
                // perform_rma = PerzynaMCCHardRMA(p_stress, q_stress, exit, M, p0, beta, xi, mu, K, dt, dim, perzyna_visc, particles.eps_pl_vol_mcc[p]);
                perform_rma = PerzynaMCCHardRMA(p_stress, q_stress, exit, M, p0, beta, xi, mu, K, dt, dim, perzyna_visc, particles.eps_pl_vol[p]);
            }
            else if (plastic_model == SinterMCC)
            {
                perform_rma = SinterMCCRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, dt, sinter_Sc, sinter_tc, sinter_ec, particles.eps_pl_vol[p], particles.sinter_S[p]);
            }

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T eps_pl_dev_instant = (q_trial - q_stress) / (mu*sqrt6);
                particles.eps_pl_dev[p] += eps_pl_dev_instant;
                particles.delta_gamma[p] = eps_pl_dev_instant / dt;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;
                particles.eps_pl_vol_mcc[p] += eps_pl_vol_inst;

                // T delta_S = dt/sinter_tc*(sinter_Sinf-particles.sinter_S[p]) - particles.sinter_S[p] * std::abs(eps_pl_vol_inst) / sinter_ec;
                // particles.sinter_S[p] = std::min(T(sinter_Sinf), std::max(T(0.0), particles.sinter_S[p] + delta_S));

                hencky = q_stress / (mu*sqrt6) * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            }
        }

    } // end plastic_model type

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

} // end plasticity
