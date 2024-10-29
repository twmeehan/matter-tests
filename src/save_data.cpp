#include "simulation.hpp"
#include "../deps/tinyply.h"

void Simulation::saveInfo(){

    std::ofstream infoFile(directory + sim_name + "/info.txt");
    infoFile << end_frame           << "\n"   // 0
             << fps                 << "\n"   // 1
             << dx                  << "\n"   // 2
             << mu                  << "\n"   // 3
             << lambda              << "\n"   // 4
             << Np                  << "\n"   // 5
             << particle_volume     << "\n"   // 6
             << rho                 << "\n"   // 7
             << grain_diameter      << "\n"   // 8
             << in_numb_ref         << "\n"   // 9
             << rho_s               << "\n"   // 10
             << mu_1                << "\n"   // 11
             << mu_2                << "\n";  // 12
    infoFile.close();
}

void Simulation::computeAvgData(TM& volavg_cauchy, TM& volavg_kirchh, T& Javg){

    volavg_cauchy = TM::Zero();
    volavg_kirchh = TM::Zero();

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

    volavg_cauchy /= Jsum;
    volavg_kirchh /= Jsum;

    Javg = Jsum / Np;

}

void Simulation::saveAvgData(){

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
    T Javg;

    computeAvgData(volavg_cauchy, volavg_kirchh, Javg);

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

    std::ofstream outFile3(directory + sim_name + "/out_Javg_frame_" + std::to_string(frame) + ".csv");
    outFile3 << Javg  << "\n";
    outFile3.close();

    std::ofstream outFile4(directory + sim_name + "/last_written.txt");
    outFile4 << std::to_string(frame) << "\n";
    outFile4.close();

}


void Simulation::saveParticleData(std::string extra){

    std::vector<T> pressure_vec(Np);
    std::vector<T> devstress_vec(Np);
    std::vector<T> Je_vec(Np);

    #pragma omp parallel for num_threads(n_threads)
    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        TM tau; // particles.tau[p];
        if (elastic_model == NeoHookean)
            tau = NeoHookeanPiola(Fe) * Fe.transpose();
        else if (elastic_model == StvkWithHencky)
            tau = StvkWithHenckyPiola(Fe) * Fe.transpose();

        T Je = Fe.determinant();

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * TM::Identity();

        T devstress;
        if (use_von_mises_q)
            devstress = std::sqrt(3.0/2.0 * selfDoubleDot(tau_dev));
        else
            devstress = std::sqrt(1.0/2.0 * selfDoubleDot(tau_dev));

        pressure_vec[p] = pressure;
        devstress_vec[p] = devstress;
        Je_vec[p] = Je;
    }

#ifdef TINYPLY_IMPLEMENTATION

    std::string filename = directory + sim_name + "/particles_f" + extra + std::to_string(frame) + ".ply";
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

#endif // TINYPLY_IMPLEMENTATION

    std::ofstream outFile(directory + sim_name + "/last_written.txt");
    outFile << std::to_string(frame) << "\n";
    outFile.close();

} // end saveParticleData()

void Simulation::saveGridData(std::string extra){

#ifdef TINYPLY_IMPLEMENTATION
    std::string filename = directory + sim_name + "/grid_f" + extra + std::to_string(frame) + ".ply";
    std::ofstream out;
    out.open(filename, std::ios::out | std::ios::binary);
    tinyply::PlyFile file;
    auto type = std::is_same<T, float>::value ? tinyply::Type::FLOAT32 : tinyply::Type::FLOAT64;

    int counter = 0;
    std::vector<TV> grid_x_save;
    grid_x_save.resize(grid_nodes);
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            #ifdef THREEDIM
            for(int k = 0; k < Nz; k++){
                TV grid_x_node(grid.x[i],grid.y[j],grid.z[k]);
            #else
                TV grid_x_node(grid.x[i],grid.y[j]);
            #endif
                grid_x_save[counter] = grid_x_node;
                counter++;
            #ifdef THREEDIM
            }
            #endif
        }
    }

    #ifdef THREEDIM
    file.add_properties_to_element(
        "vertex",
        { "x", "y", "z"},
        type,
        grid_x_save.size(),
        reinterpret_cast<uint8_t*>(grid_x_save.data()),
        tinyply::Type::INVALID,
        0);
    file.add_properties_to_element(
        "vertex",
        { "vx", "vy", "vz" },
        type,
        grid.v.size(),
        reinterpret_cast<uint8_t*>(grid.v.data()),
        tinyply::Type::INVALID,
        0);
    #else
    file.add_properties_to_element(
        "vertex",
        { "x", "y"},
        type,
        grid_x_save.size(),
        reinterpret_cast<uint8_t*>(grid_x_save.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "vx", "vy"},
        type,
        grid.v.size(),
        reinterpret_cast<uint8_t*>(grid.v.data()),
        tinyply::Type::INVALID,
        0);
    #endif

    file.add_properties_to_element(
        "vertex",
        { "mass" },
        type,
        grid.mass.size(),
        reinterpret_cast<uint8_t*>(grid.mass.data()),
        tinyply::Type::INVALID,
        0);

    if (use_material_fricton){
        file.add_properties_to_element(
            "vertex",
            { "friction" },
            type,
            grid.friction.size(),
            reinterpret_cast<uint8_t*>(grid.friction.data()),
            tinyply::Type::INVALID,
            0);
    }

    file.write(out, true);

#else

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
#endif // THREEDIM

#endif // TINYPLY_IMPLEMENTATION

} // end saveGridData()
