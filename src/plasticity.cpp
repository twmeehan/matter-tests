// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include "plasticity_helpers/mccrma.hpp"
#include "plasticity_helpers/mcchardexprma.hpp"

void Simulation::plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial){

    if (plastic_model == NoPlasticity){
        // Do nothing
    }

    else if (plastic_model == VM || plastic_model == DP || plastic_model == DPSoft || plastic_model == MCC || plastic_model == VMVisc || plastic_model == DPVisc || plastic_model == MCCVisc || plastic_model == DPMui || plastic_model == MCCMui){

        Eigen::JacobiSVD<TM> svd(Fe_trial, Eigen::ComputeFullU | Eigen::ComputeFullV);
        // TV hencky = svd.singularValues().array().log();
        TV hencky = svd.singularValues().array().abs().max(1e-4).log();
        T  hencky_trace = hencky.sum();
        TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
        T  hencky_deviatoric_norm = hencky_deviatoric.norm();

        if (hencky_deviatoric_norm > 0)
            hencky_deviatoric /= hencky_deviatoric_norm; // normalize the deviatoric vector so it gives a unit vector specifying the deviatoric direction


        if (plastic_model == VM){

            T yield_stress = q_max;

            T delta_gamma = hencky_deviatoric_norm - yield_stress / e_mu_prefac; // this is eps_pl_dev_instant

            if (delta_gamma > 0){ // project to yield surface
                plastic_count++;
                particles.delta_gamma[p] = d_prefac * delta_gamma / dt;

                hencky -= delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;

            } // end plastic projection

        } // end VonMises

        else if (plastic_model == DP){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = e_mu_prefac * hencky_deviatoric_norm;

            T q_yield = dp_slope * p_trial + dp_cohesion;

            // left of tip
            if (q_yield < 1e-10){
                plastic_count++;

                T p_proj = -dp_cohesion/dp_slope; // larger than p_trial
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();

                T delta_gamma = d_prefac * hencky_deviatoric_norm;
                particles.delta_gamma[p] = delta_gamma / dt;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.eps_pl_vol[p] += (p_proj-p_trial)/K;
            }
            else{ // right of tipe
                T delta_gamma = d_prefac * (hencky_deviatoric_norm - q_yield / e_mu_prefac);

                if (delta_gamma > 0){ // project to yield surface
                    plastic_count++;
                    particles.delta_gamma[p] = delta_gamma / dt;

                    hencky -= (1.0/d_prefac) * delta_gamma * hencky_deviatoric;
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                    particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                }
            } // if else side of tip

        } // end DP

        else if (plastic_model == DPSoft){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = e_mu_prefac * hencky_deviatoric_norm;

            T p_tip_orig = -dp_cohesion/dp_slope;
            T p_tip      = p_tip_orig * std::exp(-xi * particles.eps_pl_dev[p]);
            T p_shift = 0;
            if (use_pradhana)
                p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain!

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + (-p_tip*dp_slope); // not sure if we should really shift this intersection!!!

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) <= p_tip){
                plastic_count++;

                T p_proj = p_tip_orig * std::exp(-xi * particles.eps_pl_dev[p]); // > p_trial
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();

                T delta_gamma = d_prefac * hencky_deviatoric_norm;
                particles.delta_gamma[p] = delta_gamma / dt;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;

                T eps_pl_vol_inst = (p_proj-p_trial)/K;
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                if (use_pradhana)
                    particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
            }

            // right of tip AND inside yield surface, i.e., elastic states
            if ((p_trial+p_shift) > p_tip && q_trial <= q_yield) {
                if (use_pradhana)
                    particles.eps_pl_vol_pradhana[p] = 0; // reset correction as no longer dilating
                particles.delta_gamma[p] = 0; // elastic particles have no delta_gamma
            }

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {
                plastic_count++;

                T temp_eps_pl_dev = particles.eps_pl_dev[p] + (hencky_deviatoric_norm - q_yield / e_mu_prefac);
                T p_proj          = p_tip_orig * std::exp(-xi * temp_eps_pl_dev);

                // if left of tip: project from p_trial to p_proj
                if ((p_trial+p_shift) < p_proj){
                    T delta_gamma             = d_prefac * hencky_deviatoric_norm;
                    particles.delta_gamma[p]  = delta_gamma / dt;
                    particles.eps_pl_dev[p]  += (1.0/d_prefac) * delta_gamma;

                    T eps_pl_vol_inst = (p_proj-p_trial)/K;
                    particles.eps_pl_vol[p] += eps_pl_vol_inst;
                    if (use_pradhana)
                        particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!

                    hencky = -p_proj/(K*dim) * TV::Ones();
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                }
                // if right of tip: p_trial becomes p
                else{
                    T q_yield_new = dp_slope * (p_trial+p_shift) + (-p_proj*dp_slope) ;
                    T delta_gamma = d_prefac * (hencky_deviatoric_norm - q_yield_new / e_mu_prefac);
                    hencky -= (1.0/d_prefac) * delta_gamma * hencky_deviatoric;
                    particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                    particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                    particles.delta_gamma[p] = delta_gamma / dt;
                    if (use_pradhana)
                        particles.eps_pl_vol_pradhana[p] = 0; // reset volume accumulation
                }
            } // end plastic projection projection

        } // end DPSoft

        else if (plastic_model == VMVisc){

            // trial q-stress
            T stress = e_mu_prefac * hencky_deviatoric_norm;

            T yield_stress = q_min + (q_max - q_min) * exp(-xi * particles.eps_pl_dev[p]);

            //////////////// Only for capped von Mises /////////////////
            T p_trial = -K * hencky_trace;
            if (p_trial < p_min * exp(-xi * particles.eps_pl_vol[p])){
                T delta_gamma = stress / f_mu_prefac;
                T eps_pl_vol_inst = -p_trial/K;
                particles.F[p] = svd.matrixU() * svd.matrixV().transpose();
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;
            }
            ////////////////////////////////////////////////////////////

            else if (stress > yield_stress) {
                plastic_count++;
                T delta_gamma;

                if (perzyna_exp == 1 && xi == 0){
                    delta_gamma = (stress-yield_stress) / (f_mu_prefac + yield_stress*perzyna_visc/dt);
                    if (delta_gamma < 0){
                        debug("VMVisc: FATAL negative delta_gamma = ", delta_gamma);
                        exit = 1;
                    }
                }
                else{

                    delta_gamma = 0.1 * (stress - yield_stress) / f_mu_prefac; // initial guess
                    int max_iter = 60;
                    for (int iter = 0; iter < max_iter; iter++) {
                        if (iter == max_iter - 1){ // did not break loop
                            debug("VMVisc: FATAL did not exit loop at iter = ", iter);
                            exit = 1;
                        }

                        if (delta_gamma < 0) // not possible and can also lead to division by zero
                            delta_gamma = 1e-10;

                        T tm = perzyna_visc * delta_gamma + dt;
                        T tmp = dt / tm;
                        T tmp1 = std::pow(tmp, perzyna_exp);

                        T yield_stress_new = q_min + (q_max - q_min) * exp(-xi * (particles.eps_pl_dev[p] + (1.0/d_prefac) * delta_gamma));

                        T residual = (stress - f_mu_prefac * delta_gamma) * tmp1 - yield_stress_new;
                        if (std::abs(residual) < 1e-1) {
                            break;
                        }

                        T yield_stress_new_diff = -xi / d_prefac * (q_max - q_min) * exp(-xi * (particles.eps_pl_dev[p] + (1.0/d_prefac) * delta_gamma));
                        T residual_diff         = -f_mu_prefac * tmp1 + (stress - f_mu_prefac * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm) - yield_stress_new_diff;

                        if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                            debug("VMVisc: residual_diff too small in abs value = ", residual_diff);
                            exit = 1;
                        }

                        delta_gamma -= residual / residual_diff;
                    } // end N-R iterations
                } // end if perzyna_exp > 1

                hencky -= (1.0/d_prefac) * delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;
            } // end plastic projection projection
        } // end VMVisc

        else if (plastic_model == DPVisc){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = e_mu_prefac * hencky_deviatoric_norm;

            T p_tip   = -dp_cohesion/dp_slope;
            T p_shift = 0;
            if (use_pradhana)
                p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            // if positive volume gain, the q=0 intersection for the yield surface is shifted to the right, at a larger p.
            T q_yield = dp_slope * (p_trial+p_shift) + dp_cohesion;

            if (use_material_friction)
                particles.muI[p] = dp_slope;

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) < p_tip){
                T delta_gamma     = d_prefac * hencky_deviatoric_norm;
                T p_proj          = p_tip; // > p_trial
                T eps_pl_vol_inst = (p_proj-p_trial)/K; // this can be both positive and negative when using Pradhana!
                plastic_count++;
                hencky = -p_proj/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.delta_gamma[p] = delta_gamma / dt;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.eps_pl_vol[p] += eps_pl_vol_inst;
                if (use_pradhana)
                    particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
            }
            else{ // if right of shifted tip (incl elastic states)
                if (use_pradhana)
                    particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
                particles.delta_gamma[p] = 0; // for the elastic particles, the plastic particles have their delta_gamma overwritten in the next if-statement
            }

            // right of shifted tip AND outside the shifted yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > q_yield) {
                plastic_count++;

                T delta_gamma;

                if (perzyna_exp == 1){
                    delta_gamma = (q_trial-q_yield) / (f_mu_prefac + q_yield*perzyna_visc/dt);
                    if (delta_gamma < 0){
                        debug("DPVisc: FATAL negative delta_gamma = ", delta_gamma);
                        exit = 1;
                    }
                }
                else{

                    delta_gamma = 0.01 * (q_trial - q_yield) / f_mu_prefac; // initial guess

                    int max_iter = 60;
                    for (int iter = 0; iter < max_iter; iter++) {
                        if (iter == max_iter - 1){ // did not break loop
                            debug("DPVisc: FATAL did not exit loop at iter = ", iter);
                            exit = 1;
                        }

                        T tm = perzyna_visc * delta_gamma + dt;
                        T tmp = dt / tm;
                        T tmp1 = std::pow(tmp, perzyna_exp);

                        T residual = (q_trial - f_mu_prefac * delta_gamma) * tmp1 - q_yield;
                        if (std::abs(residual) < 1e-2) {
                            break;
                        }

                        T residual_diff = -f_mu_prefac * tmp1 + (q_trial - f_mu_prefac * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm);

                        if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                            debug("DPVisc: FATAL residual_diff too small in abs value = ", residual_diff);
                            exit = 1;
                        }

                        delta_gamma -= residual / residual_diff;

                        if (delta_gamma < 0) // not possible and can also lead to division by zero
                            delta_gamma = 1e-10;

                    } // end N-R iterations

                } // end if perzyna_exp == 1

                if (use_material_friction){
                    if (std::abs(p_trial + p_shift) < 1e-10)
                        particles.muI[p] = 1e15;
                    else
                        particles.muI[p] = ((q_trial - f_mu_prefac * delta_gamma) - dp_cohesion) / (p_trial + p_shift);
                }

                hencky -= (1.0/d_prefac) * delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;
            } // end plastic projection
        } // end DPVisc


        else if (plastic_model == DPMui){

            // trial stresses
            T p_trial = -K * hencky_trace;
            T q_trial = e_mu_prefac * hencky_deviatoric_norm;

            T p_tip   = -dp_cohesion/mu_1;
            T p_shift = 0;
            if (use_pradhana)
                p_shift = -K * particles.eps_pl_vol_pradhana[p]; // Negative if volume gain! Force to be zero if using classical volume-expanding non-ass. DP

            particles.muI[p]         = mu_1;
            particles.viscosity[p]   = 0;
            particles.delta_gamma[p] = 0;

            // if left of shifted tip,
            // => project to the original tip given by cohesion only (i.e., not the shifted tip)
            if ((p_trial+p_shift) < p_tip){
                T delta_gamma     = d_prefac * hencky_deviatoric_norm;
                T eps_pl_vol_inst = (p_tip-p_trial)/K;
                plastic_count++;
                hencky = -p_tip/(K*dim) * TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.delta_gamma[p]          = delta_gamma;               // NB!
                particles.eps_pl_dev[p]          += (1.0/d_prefac) * delta_gamma;
                particles.eps_pl_vol[p]          += eps_pl_vol_inst;
                particles.eps_pl_vol_pradhana[p] += eps_pl_vol_inst; // can be negative!
            }
            else{ // if right of shifted tip (incl elastic states)
                particles.eps_pl_vol_pradhana[p] = 0; // reset pradhana volume accumulation
            }

            // if positive volume gain, the q=0 intersection for the plastic potential surface is shifted to the right, at a larger p.
            T q_yield = mu_1 * (p_trial+p_shift) + dp_cohesion;

            // right of tip AND outside yield surface
            if ((p_trial+p_shift) > p_tip && q_trial > (q_yield + 1e-4)) {

                plastic_count++;

                T p_special = p_trial+p_shift - p_tip; // p_special is positive.

                T fac_a = f_mu_prefac * dt; // always positive
                T fac_b = p_trial*(mu_2-mu_1) + f_mu_prefac*dt*fac_Q*std::sqrt(p_special) - (q_trial-q_yield);
                T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(p_special); // always negative (or zero IFF p_special = 0)

                // this is gamma_dot:
                T delta_gamma = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always if because a>0 and c<0

                T mu_i = mu_2; // if p_special = 0 then mu = mu_2
                if (p_special > 1e-4){
                    mu_i = mu_1 + (mu_2 - mu_1) / (fac_Q * std::sqrt(p_special) / delta_gamma + 1.0);
                }
                particles.muI[p]       = mu_i;
                particles.viscosity[p] = (mu_i - mu_1) * p_special / delta_gamma;

                delta_gamma *= dt; // this is the actual delta_gamma

                hencky -= (1.0/d_prefac) * delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;

            } // end plastic projection

        } // end DPMui

        else if (plastic_model == MCCMui) {

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = e_mu_prefac * hencky_deviatoric_norm;

            // make copies
            T q_trial = q_stress;

            //////////////////////////////////////////////////////////////////////
            bool perform_rma;
            T p_c;
            if (hardening_law == ExpExpl){ // Exponential Explicit Hardening
                p_c = std::max(T(1e-2), p0*std::exp(-xi*particles.eps_pl_vol[p]));
                perform_rma = MCCRMA(p_stress, q_stress, exit, mu_1, p_c, beta, mu, K, rma_prefac);
            }
            else if (hardening_law == SinhExpl){ // Sinh Explicit Hardening
                p_c = std::max(T(1e-2), K*std::sinh(-xi*particles.eps_pl_vol[p] + std::asinh(p0/K)));
                perform_rma = MCCRMA(p_stress, q_stress, exit, mu_1, p_c, beta, mu, K, rma_prefac);
            }
            else if (hardening_law == ExpImpl){ // Exponential Implicit Hardening
                perform_rma = MCCHardExpRMA(p_stress, q_stress, exit, mu_1, p0, beta, mu, K, xi, rma_prefac, particles.eps_pl_vol[p]);
            }
            else{
                debug("You specified an unvalid HARDENING LAW!");
                exit = 1;
            }

            //////////////////////////////////////////////////////////////////////

            particles.muI[p] = mu_1;
            particles.viscosity[p] = 0;

            if (perform_rma){
                plastic_count++;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p]     += eps_pl_vol_inst;

                T q_yield = q_stress;
                if (q_trial < (q_yield + 1e-4)) {
                    q_stress = q_yield; // this is to ensure that fac_c below becomes negative and thus gamma_dot_S positive
                }
                else{
                    if (hardening_law == ExpImpl)
                        p_c = std::max( p0 * std::exp(-xi*particles.eps_pl_vol[p]) , p_stress+T(0.001) );

                    T p_special = p_stress + beta * p_c; // always equal or larger than 0

                    T fac_a = f_mu_prefac * dt; // always positive

                    T fac_b = std::sqrt((p_c-p_stress)*p_special)*(mu_2-mu_1) + f_mu_prefac*dt*fac_Q*std::sqrt(p_special) - (q_trial-q_yield);

                    T fac_c = -(q_trial-q_yield) * fac_Q * std::sqrt(p_special); // always negative

                    T gamma_dot_S = (-fac_b + std::sqrt(fac_b*fac_b - 4*fac_a*fac_c) ) / (2*fac_a); // always positive because a>0 and c<0

                    q_stress = std::max(q_yield, q_trial - f_mu_prefac * dt * gamma_dot_S);

                    T mu_i = mu_2; // if p_special = 0 then mu = mu_2
                    if (p_special > 1e-5){
                        mu_i = mu_1 + (mu_2 - mu_1) / (fac_Q * std::sqrt(p_special) / gamma_dot_S + 1.0);
                    }
                    particles.muI[p]       = mu_i;
                    particles.viscosity[p] = (mu_i - mu_1) * std::sqrt((p_c-p_stress)*p_special) / gamma_dot_S;
                }

                T delta_gamma = (q_trial - q_stress) / f_mu_prefac;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;

                hencky = q_stress / e_mu_prefac * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            } // end plastic projection

        } // end MCCMui


        else if (plastic_model == MCC || plastic_model == MCCVisc){

            // the trial stress states
            T p_stress = -K * hencky_trace;
            T q_stress = e_mu_prefac * hencky_deviatoric_norm;

            // make copies
            T p_trial = p_stress;
            T q_trial = q_stress;

            bool perform_rma;

            if (hardening_law == NoHard){ // Exponential Explicit Hardening
                perform_rma = MCCRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, rma_prefac);
            }
            else if (hardening_law == ExpExpl){ // Exponential Explicit Hardening
                T particle_p0 = std::max(T(1e-2), p0*std::exp(-xi*particles.eps_pl_vol[p]));
                // T particle_p0 = std::max( T(1e-2), p0*std::exp(xi*(1-std::exp(particles.eps_pl_vol[p]))) );
                perform_rma = MCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, rma_prefac);
            }
            else if (hardening_law == SinhExpl){ // Sinh Explicit Hardening
                T particle_p0 = std::max(T(1e-2), K*std::sinh(-xi*particles.eps_pl_vol[p] + std::asinh(p0/K)));
                // T particle_p0 = std::max(T(1e-2), p0*std::sinh(-xi*particles.eps_pl_vol[p] + std::asinh(1.0)));
                // T particle_p0 = std::max(T(1e-2), (particles.eps_pl_vol[p] < 0) ? p0*(1.0-std::sinh(xi*particles.eps_pl_vol[p])) : p0*(1.0-std::tanh(xi*particles.eps_pl_vol[p])) );
                perform_rma = MCCRMA(p_stress, q_stress, exit, M, particle_p0, beta, mu, K, rma_prefac);
            }
            else if (hardening_law == ExpImpl){ // Exponential Implicit Hardening
                perform_rma = MCCHardExpRMA(p_stress, q_stress, exit, M, p0, beta, mu, K, xi, rma_prefac, particles.eps_pl_vol[p]);
            }
            else{
                debug("You specified an unvalid HARDENING LAW!");
                exit = 1;
            }

            if (perform_rma) { // returns true if it performs a return mapping
                plastic_count++;

                T ep = p_stress / (K*dim);
                T eps_pl_vol_inst = hencky_trace + dim * ep;
                particles.eps_pl_vol[p] += eps_pl_vol_inst;

                T delta_gamma;
                ////////////////////////////////////////////////////////////////
                if (plastic_model == MCCVisc){

                    T q_yield = q_stress;

                    if (perzyna_exp == 1){
                        delta_gamma = (q_trial-q_yield) / (f_mu_prefac + q_yield*perzyna_visc/dt);
                        if (delta_gamma < 0){
                            debug("MCCVisc: FATAL negative delta_gamma = ", delta_gamma);
                            // delta_gamma = 0;
                            exit = 1;
                        }
                    }
                    else{ // persyna_exp is not one

                        delta_gamma = 0.01 * (q_trial - q_yield) / f_mu_prefac; // initial guess

                        int max_iter = 60;
                        for (int iter = 0; iter < max_iter; iter++) {
                            if (iter == max_iter - 1){ // did not break loop
                                debug("MCCVisc: FATAL did not exit loop at iter = ", iter);
                                exit = 1;
                            }

                            T tm = perzyna_visc * delta_gamma + dt;
                            T tmp = dt / tm;
                            T tmp1 = std::pow(tmp, perzyna_exp);

                            T residual = (q_trial - f_mu_prefac * delta_gamma) * tmp1 - q_yield;
                            if (std::abs(residual) < 1e-2) {
                                break;
                            }

                            T residual_diff = -f_mu_prefac * tmp1 + (q_trial - f_mu_prefac * delta_gamma) * perzyna_exp * std::pow(tmp, perzyna_exp - 1) * (-perzyna_visc * dt) / (tm * tm);

                            if (std::abs(residual_diff) < 1e-14){ // otherwise division by zero
                                debug("MCCVisc: FATAL residual_diff too small in abs value = ", residual_diff);
                                exit = 1;
                            }

                            delta_gamma -= residual / residual_diff;

                            if (delta_gamma < 0) // not possible and can also lead to division by zero
                                delta_gamma = 1e-10;

                        } // end N-R iterations

                    } // end if perzyna_exp == 1

                    q_stress = std::max(q_yield, q_trial - f_mu_prefac * delta_gamma); // delta_gamma = dt * gamma_dot_S

                } //  end if MCCVisc
                ////////////////////////////////////////////////////////////////


                delta_gamma = (q_trial - q_stress) / f_mu_prefac;
                particles.eps_pl_dev[p] += (1.0/d_prefac) * delta_gamma;
                particles.delta_gamma[p] = delta_gamma / dt;

                hencky = q_stress / e_mu_prefac * hencky_deviatoric - ep*TV::Ones();
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            } // if perform_rma
        } // end MCC / MCCVisc

    } // end plastic_model type

    else{
        debug("You specified an unvalid PLASTIC model!");
    }

} // end plasticity
