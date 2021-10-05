#include "simulation.hpp"

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
            << "delta_gamma"         << ","   // 10
            << "eps_pl_vol_pradhana" << "\n";   // 11
            // << "tau_xx"      << ","   // 13
            // << "tau_xy"      << ","   // 14
            // << "tau_yx"      << ","   // 15
            // << "tau_yy"      << ","   // 16
            // << "Fe_xx"       << ","   // 17
            // << "Fe_xy"       << ","   // 18
            // << "Fe_yx"       << ","   // 19
            // << "Fe_yy"       << "\n";   // 20

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

        T J = Fe.determinant() * std::exp( particles.eps_pl_vol[p] );
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
                << particles.delta_gamma[p]          << ","     // 10
                << particles.eps_pl_vol_pradhana[p]  << "\n";   // 11
                // << tau(0,0)                   << ","   // 13
                // << tau(0,1)                   << ","   // 14
                // << tau(1,0)                   << ","   // 15
                // << tau(1,1)                   << ","   // 16
                // << Fe(0,0)                    << ","    // 17
                // << Fe(0,1)                    << ","    // 18
                // << Fe(1,0)                    << ","    // 19
                // << Fe(1,1)                    << "\n";  // 20
    } // end loop over particles
    outFile.close();

    volavg_tau /= Jsum;
    T volavg_pressure = -volavg_tau.trace() / dim;
    TM volavg_tau_dev = volavg_tau + volavg_pressure * I;
    T volavg_devstress = std::sqrt(3.0/2.0 * selfDoubleDot(volavg_tau_dev));
    std::ofstream outFile2(directory + sim_name + "/out_pq_frame_" + extra + std::to_string(frame) + ".csv");
    outFile2 << volavg_pressure    << ","
             << volavg_devstress   << "\n";
    outFile2.close();

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
