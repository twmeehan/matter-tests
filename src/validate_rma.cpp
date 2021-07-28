#include "simulation.hpp"

void Simulation::validateRMA(){

    createDirectory();

    std::ofstream plastic_info; plastic_info.open(directory + sim_name + "/plastic_info.txt");
    plastic_info << M << "\n" << p0 << "\n" << beta << std::endl;
    plastic_info.close();

    // Trial state ENTER HERE
    T p_trial = 0.345086;
    T q_trial = 0.720833;

    T trace_epsilon = -p_trial / K;
    T norm_eps_hat = q_trial / mu_sqrt6;

    std::ofstream steps; steps.open(directory + sim_name + "/rma_steps.txt");
    steps    << "0" << "\t" << p_trial << "\t" << q_trial << "\t" << "0" << std::endl;

    bool outside = AnalQuadReturnMapping(p_trial, q_trial, exit, M, p0, beta);
    // bool outside = QuadraticReturnMapping(p_trial, q_trial, exit, trace_epsilon, norm_eps_hat, M, p0, beta, mu, K);
    // bool outside = CamClayReturnMapping(p_trial, q_trial, exit, trace_epsilon, norm_eps_hat, M, p0, beta, mu, K);

    if (exit == 1){
        std::cout << "RMA failed with exit  = " << exit << std::endl;
        return;
    }

    steps    << "1" << "\t" << p_trial << "\t" << q_trial << "\t" << "0" << std::endl;
    steps.close();

    std::cout << "outside  = " << outside << std::endl;
    std::cout << "p_proj  = " << p_trial << std::endl;
    std::cout << "q_proj  = " << q_trial << std::endl;

    std::cout << "For this p, the correct q is " << 2.0*M / (2*beta+1.0) * (p0-p_trial)*(beta*p0+p_trial) / p0 << std::endl;
    //std::cout << "For this p, the correct q is " << M * std::sqrt( (p0-p_trial)*(sbeta*p0+p_trial) / (2*beta+1.0) ) << std::endl;

    std::cout << "---------------------------"  << std::endl;

}
