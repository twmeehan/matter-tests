#include "simulation.hpp"
#include "tinyply.h"

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
             << rho                 << "\n"   // 9
             << grain_diameter      << "\n"   // 10
             << in_numb_ref         << "\n"   // 11
             << rho_s               << "\n"   // 12
             << mu_1                << "\n"   // 13
             << mu_2                << "\n";  // 14
    infoFile.close();
}

void Simulation::saveAvgData(){

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
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
        Jsum += J;

        volavg_cauchy += tau;
        volavg_kirchh += tau * J;
    }

    // T volavg_pressure = -volavg_tau.trace() / dim;
    // TM volavg_tau_dev = volavg_tau + volavg_pressure * I;
    // T volavg_devstress = std::sqrt(3.0/2.0 * selfDoubleDot(volavg_tau_dev));

    volavg_cauchy /= Np;
    volavg_kirchh /= Np;
    // must later be multiplied by phi_0/(1+eps_V)

    std::ofstream outFile1(directory + sim_name + "/out_avgcauchy_frame_" + std::to_string(frame) + ".csv");
    outFile1 << volavg_cauchy(0,0)    << ","
             << volavg_cauchy(0,1)    << ","
             << volavg_cauchy(1,0)    << ","
             << volavg_cauchy(1,1)    << "\n";
    outFile1.close();

    std::ofstream outFile2(directory + sim_name + "/out_avgkirchh_frame_" + std::to_string(frame) + ".csv");
    outFile2 << volavg_kirchh(0,0)    << ","
             << volavg_kirchh(0,1)    << ","
             << volavg_kirchh(1,0)    << ","
             << volavg_kirchh(1,1)    << "\n";
    outFile2.close();

    std::ofstream outFile3(directory + sim_name + "/out_avgdetF_frame_" + std::to_string(frame) + ".csv");
    outFile3 << Jsum/Np  << "\n";
    outFile3.close();

    std::ofstream outFile4(directory + sim_name + "/last_written.txt");
    outFile4 << std::to_string(frame) << "\n";
    outFile4.close();

}


void Simulation::saveParticleData(std::string extra){

    std::vector<T> pressure_vec(Np);
    std::vector<T> devstress_vec(Np);
    std::vector<T> Je_vec(Np);

    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        TM tau; // particles.tau[p];
        if (elastic_model == NeoHookean)
            tau = NeoHookeanPiola(Fe) * Fe.transpose();
        else if (elastic_model == StvkWithHencky)
            tau = StvkWithHenckyPiola(Fe) * Fe.transpose();

        T Je = Fe.determinant();
        // T J = Je * std::exp( particles.eps_pl_vol[p] );

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * TM::Identity();
        T devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));

        pressure_vec[p] = pressure;
        devstress_vec[p] = devstress;
        Je_vec[p] = Je;
    }

#ifdef TINYPLY_IMPLEMENTATION

    std::string filename = directory + sim_name + "/out_part_frame_" + extra + std::to_string(frame) + ".ply";
    std::ofstream out;
    out.open(filename, std::ios::out | std::ios::binary);
    tinyply::PlyFile file;

    auto type = std::is_same<T, float>::value ? tinyply::Type::FLOAT32 : tinyply::Type::FLOAT64;

    #ifdef THREEDIM
    file.add_properties_to_element(
        "vertex",
        { "x", "y", "z"},
        type,
        particles.x.size(),
        reinterpret_cast<uint8_t*>(particles.x.data()),
        tinyply::Type::INVALID,
        0);
    file.add_properties_to_element(
        "vertex",
        { "vx", "vy", "vz" },
        type,
        particles.v.size(),
        reinterpret_cast<uint8_t*>(particles.v.data()),
        tinyply::Type::INVALID,
        0);
    #else
    file.add_properties_to_element(
        "vertex",
        { "x", "y"},
        type,
        particles.x.size(),
        reinterpret_cast<uint8_t*>(particles.x.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "vx", "vy"},
        type,
        particles.v.size(),
        reinterpret_cast<uint8_t*>(particles.v.data()),
        tinyply::Type::INVALID,
        0);
    #endif

    file.add_properties_to_element(
        "vertex",
        { "eps_pl_vol" },
        type,
        particles.eps_pl_vol.size(),
        reinterpret_cast<uint8_t*>(particles.eps_pl_vol.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "eps_pl_dev" },
        type,
        particles.eps_pl_dev.size(),
        reinterpret_cast<uint8_t*>(particles.eps_pl_dev.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "delta_gamma" },
        type,
        particles.delta_gamma.size(),
        reinterpret_cast<uint8_t*>(particles.delta_gamma.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "viscosity" },
        type,
        particles.viscosity.size(),
        reinterpret_cast<uint8_t*>(particles.viscosity.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "muI" },
        type,
        particles.muI.size(),
        reinterpret_cast<uint8_t*>(particles.muI.data()),
        tinyply::Type::INVALID,
        0);

    // file.add_properties_to_element(
    //     "vertex",
    //     { "eps_pl_vol_pradhana" },
    //     type,
    //     particles.eps_pl_vol_pradhana.size(),
    //     reinterpret_cast<uint8_t*>(particles.eps_pl_vol_pradhana.data()),
    //     tinyply::Type::INVALID,
    //     0);
    //
    // file.add_properties_to_element(
    //     "vertex",
    //     { "sinter_S" },
    //     type,
    //     particles.sinter_S.size(),
    //     reinterpret_cast<uint8_t*>(particles.sinter_S.data()),
    //     tinyply::Type::INVALID,
    //     0);

    file.add_properties_to_element(
        "vertex",
        { "pressure" },
        type,
        pressure_vec.size(),
        reinterpret_cast<uint8_t*>(pressure_vec.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "devstress" },
        type,
        devstress_vec.size(),
        reinterpret_cast<uint8_t*>(devstress_vec.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "Je" },
        type,
        Je_vec.size(),
        reinterpret_cast<uint8_t*>(Je_vec.data()),
        tinyply::Type::INVALID,
        0);

    file.write(out, true);

#else

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
                << pressure_vec[p]            << ","   // 6
                << devstress_vec[p]           << ","   // 7
                << particles.eps_pl_vol[p]    << ","   // 8
                << particles.eps_pl_dev[p]    << ","   // 9
                << particles.delta_gamma[p]   << ","     // 10
                << particles.viscosity[p]     << ","     // 11
                << particles.muI[p]           << ","     // 12
                << particles.eps_pl_vol_pradhana[p]  << ","  // 13
                << Je_vec[p]                 << ","
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

#endif

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
