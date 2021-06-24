#include "simulation.hpp"

void Simulation::plasticity_projection(){

    P2G_nonlocal();
    G2P_nonlocal();

    int plastic_count = 0;

    #pragma omp parallel for num_threads(n_threads)
    for(int p=0; p<Np; p++){

        T delta_gamma = particles.delta_gamma[p];

        if (delta_gamma > 0){ // project to yield surface
            plastic_count++;

            TV hencky = particles.hencky[p];
            T  hencky_trace = hencky.sum();
            TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
            T  hencky_deviatoric_norm = hencky_deviatoric.norm();

            hencky -= delta_gamma * (hencky_deviatoric / hencky_deviatoric_norm);
            Eigen::JacobiSVD<TM> svd(particles.F[p], Eigen::ComputeFullU | Eigen::ComputeFullV);
            particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
            particles.eps_pl_dev_inst[p]  = delta_gamma;
            particles.eps_pl_dev[p]      += delta_gamma;
        }

    } // end for loop

    debug("               projected nonlocal particles = ", plastic_count, " / ", Np);

} // end plasticity_projection
