#include "simulation.hpp"

void Simulation::plasticity_projection(){

    P2G_nonlocal();
    G2P_nonlocal();

    int plastic_count = 0;

    #pragma omp parallel for reduction(+:plastic_count) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        // if local approach
        // T delta_gamma_nonloc = particles.delta_gamma[p];

        // if nonlocal approach
        T delta_gamma_nonloc = particles.delta_gamma_nonloc[p];

        if (delta_gamma_nonloc > 0){ // If delta_gamma_nonloc is > 0, then delta_gamma > 0 also
            plastic_count++;

            T delta_gamma = particles.delta_gamma[p];

            if (delta_gamma > 0){ // if originally PLASTIC particle, we perform the projection using the NONLOCAL delta gamma

                TV hencky = particles.hencky[p];
                T  hencky_trace = hencky.sum();
                TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
                T  hencky_deviatoric_norm = hencky_deviatoric.norm();

                hencky -= delta_gamma_nonloc * (hencky_deviatoric / hencky_deviatoric_norm);
                Eigen::JacobiSVD<TM> svd(particles.F[p], Eigen::ComputeFullU | Eigen::ComputeFullV);
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();

                particles.eps_pl_dev[p] += delta_gamma;

            }
            // In either case of an elastic or plastic particle, we record the plastic strain. This will change the yield stress in the next time step.
            particles.eps_pl_dev_nonloc[p] += delta_gamma_nonloc;
        }

    } // end for loop

    debug("               projected nonlocal particles = ", plastic_count, " / ", Np);

} // end plasticity_projection
