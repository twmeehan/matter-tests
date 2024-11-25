// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::validateRMA(){

    createDirectory();

    std::ofstream plastic_info; plastic_info.open(directory + sim_name + "/plastic_info.txt");
    plastic_info << M << "\n" << p0 << "\n" << beta << std::endl;
    plastic_info.close();

    // Trial state ENTER HERE
    T p_trial = -100;
    T q_trial = 200;

    T trace_epsilon = -p_trial / K;
    T norm_eps_hat = q_trial / (mu*sqrt6);

    std::ofstream steps; steps.open(directory + sim_name + "/rma_steps.txt");
    steps    << "0" << "\t" << p_trial << "\t" << q_trial << "\t" << "0" << std::endl;

=    bool outside = MCCRMA(p_trial, q_trial, exit, M, p0, beta, mu, K, 1);


    if (exit == 1){
        std::cout << "RMA failed with exit  = " << exit << std::endl;
        return;
    }

    steps    << "1" << "\t" << p_trial << "\t" << q_trial << "\t" << "0" << std::endl;
    steps.close();

    std::cout << "outside  = " << outside << std::endl;
    std::cout << "p_proj  = " << p_trial << std::endl;
    std::cout << "q_proj  = " << q_trial << std::endl;

    std::cout << "For this p, the correct q is " << M * std::sqrt( (p0-p_trial)*(beta*p0+p_trial) / (2*beta+1.0) ) << std::endl;

    std::cout << "---------------------------"  << std::endl;

}
