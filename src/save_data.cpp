#include "simulation.hpp"

void Simulation::saveInfo(){

    std::ofstream infoFile(directory + sim_name + "/info.txt");
    infoFile << end_frame           << "\n"   // 0
             << fps                 << "\n"   // 1
             << dx                  << "\n"   // 2
             << mu                  << "\n"   // 3
             << lambda              << "\n"   // 4
             << vmin_factor         << "\n"   // 5
             << load_factor         << "\n"   // 6
             << Np                  << "\n"   // 7
             << particle_volume     << "\n"   // 8
             << rho                 << "\n";  // 9
    infoFile.close();
}

void Simulation::saveParticleData(std::string extra){
    std::ofstream outFile(directory + sim_name + "/out_part_frame_" + extra + std::to_string(frame) + ".csv");

    outFile << "x"           << ","   // 0
            << "y"           << ","   // 1
            << "z"           << ","   // 2
            << "vx"          << ","   // 3
            << "vy"          << ","   // 4
            << "vz"          << ","   // 5
            << "pressure"    << ","   // 6
            << "devstress"   << ","   // 7
            << "eps_pl_vol"  << ","   // 8
            << "eps_pl_dev"  << ","   // 9
            << "delta_gamma" << ","   // 10
            << "viscosity"   << ","   // 11
            << "muI"         << ","   // 12
            << "eps_pl_vol_pradhana" << ","   // 13
            << "Je"         << ","    // 14
            << "sinter_S"   << "\n";  // 15
            // << "tau_xx"      << ","   // 16
            // << "tau_xy"      << ","   // 17
            // << "tau_yx"      << ","   // 18
            // << "tau_yy"      << ","   // 19
            // << "Fe_xx"       << ","   // 20
            // << "Fe_xy"       << ","   // 21
            // << "Fe_yx"       << ","   // 22
            // << "Fe_yy"       << "\n";   // 23

    TM I = TM::Identity();
    TM volavg_tau = TM::Zero();
    T Jsum = 0;
    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        TM tau; // particles.tau[p];
        if (elastic_model == NeoHookean)
            tau = NeoHookeanPiola(Fe) * Fe.transpose();
        else if (elastic_model == StvkWithHencky)
            tau = StvkWithHenckyPiola(Fe) * Fe.transpose();

        T Je = Fe.determinant();
        T J = Je * std::exp( particles.eps_pl_vol[p] );
        volavg_tau += tau * J;
        Jsum += J;

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * I;
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        outFile << particles.x[p](0)          << ","   // 0
                << particles.x[p](1)          << ","   // 1
            #ifdef THREEDIM
                << particles.x[p](2)          << ","   // 2
            #else
                << 0                          << ","
            #endif
                << particles.v[p](0)          << ","   // 3
                << particles.v[p](1)          << ","   // 4
            #ifdef THREEDIM
                << particles.v[p](2)          << ","   // 5
            #else
                << 0                          << ","
            #endif
                << pressure                   << ","   // 6
                << devstress                  << ","   // 7
                << particles.eps_pl_vol[p]    << ","   // 8
                << particles.eps_pl_dev[p]    << ","   // 9
                << particles.delta_gamma[p]   << ","     // 10
                << particles.viscosity[p]     << ","     // 11
                << particles.muI[p]           << ","     // 12
                << particles.eps_pl_vol_pradhana[p]  << ","  // 13
                << Je                         << ","
                << particles.sinter_S[p]     << "\n";
                // << tau(0,0)                   << ","   // 14
                // << tau(0,1)                   << ","   // 15
                // << tau(1,0)                   << ","   // 16
                // << tau(1,1)                   << ","   // 17
                // << Fe(0,0)                    << ","    // 18
                // << Fe(0,1)                    << ","    // 19
                // << Fe(1,0)                    << ","    // 20
                // << Fe(1,1)                    << "\n";  // 21
    } // end loop over particles
    outFile.close();

    // UNCOMMENT IF OUTPUTTING VOL-AVG STRESS INVARIANTS:
/*
    volavg_tau /= Jsum;
    T volavg_pressure = -volavg_tau.trace() / dim;
    TM volavg_tau_dev = volavg_tau + volavg_pressure * I;
    T volavg_devstress = std::sqrt(3.0/2.0 * selfDoubleDot(volavg_tau_dev));
    std::ofstream outFile2(directory + sim_name + "/out_pq_frame_" + extra + std::to_string(frame) + ".csv");
    outFile2 << volavg_pressure    << ","
             << volavg_devstress   << ","
             << Jsum               << "\n";
    outFile2.close();
*/
    std::ofstream outFile3(directory + sim_name + "/last_written.txt");
    outFile3 << std::to_string(frame) << "\n";
    outFile3.close();

} // end saveParticleData()

void Simulation::saveGridData(std::string extra){
    std::ofstream outFile(directory + sim_name + "/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    outFile         << "x"       << ","
                    << "y"       << ","
                    << "z"       << ","
                    << "vx"      << ","
                    << "vy"      << ","
                    << "vz"      << ","
                    << "mass"    << "\n";
                //    << "delta_gamma" << ","; // OBS: Only if nonlocal strategy

#ifdef THREEDIM
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; k++){
                unsigned int index = ind(i,j,k);
                outFile << grid.x[i]               << "," // 0
                        << grid.y[j]               << "," // 1
                        << grid.z[k]               << "," // 2
                        << grid.v[index](0)        << "," // 3
                        << grid.v[index](1)        << "," // 4
                        << grid.v[index](2)        << "," // 5
                        << grid.mass[index]        << "\n";   // 6
                    //    << grid.delta_gamma[index] << ",";
            } // end for k
        } // end for j
    } // end for i
    outFile.close();
#else
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            unsigned int index = ind(i,j);
            outFile << grid.x[i]               << "," // 0
                    << grid.y[j]               << "," // 1
                    << 0                       << ","
                    << grid.v[index](0)        << "," // 3
                    << grid.v[index](1)        << "," // 4
                    << 0                       << ","
                    << grid.mass[index]        << "\n";  // 6
                //    << grid.delta_gamma[index] << ","; // 7
        } // end for j
    } // end for i
    outFile.close();
#endif

} // end saveGridData()
