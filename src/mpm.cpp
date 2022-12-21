#include "simulation.hpp"
#include "tools.hpp"
#include "sampling_particles.hpp"


// TODO:
//  * CUDA Parallilize

// REMEMBER TO SPECIFY DIMENSION AND SPLINE DEGREE IN TOOLS!!!!!

int main(){

    Simulation sim;

    sim.directory = "/media/blatny/harddrive4/larsie/";
    sim.end_frame = 100;

    unsigned int setup = 3;

    if (setup == 0){      // pbc
        sim.fps = 20;
        sim.sim_name = "pbc_kam_a28_d2-5mm_initv";
    }
    else if (setup == 1){ // gran coll
        sim.fps = 238.2166843;
        sim.sim_name = "grancoll_kam_mccmui_vgate045_new_dx2";
        // sim.sim_name = "grancoll3D_kam_mccmui_vgate045_new";
    }
    else if (setup == 2){ // incl slope
        sim.fps = 50;
        sim.sim_name = "inclslope10_mu245_kam_mccmui_comp3d";
        // sim.sim_name = "incl3Dslope22_mu245_kam_mccmui_new2";
    }
    else {                 // feeder
        sim.fps = 7;
        sim.sim_name = "feeder_kam_a24_d2-5mm_dx35_Lx18";
    }

    T theta_deg = 24;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);
    sim.gravity_time = 0.0;

    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;
    sim.n_threads = 4;

    sim.initialize(/* E */ 1e6, /* nu */ 0.3, /* rho */ 1550); // feeder, incl slope and gran collapse

    T h_gate, l_gate;
    if (setup == 0){
        sim.pbc = true;
        sim.Lx = 12*0.00125;
        sim.Ly = 0.05;
        #ifndef THREEDIM
            SampleParticles(sim.Lx, sim.Ly, 0.00023, 20, 0, sim);
        #endif
        sim.dx = 0.00125;
        sim.particle_volume = sim.Lx * sim.Ly / sim.Np;
        sim.particle_mass = sim.rho * sim.particle_volume;

        for(int p = 0; p < sim.Np; p++){
            T y = sim.particles.x[p](1);
            sim.particles.v[p](0) = y / sim.Ly * 4;
        }

        // sim.pbc = true;
        // sim.Np = 20;
        // sim.Ly = 1;
        // sim.dx = sim.Ly / sim.Np; // one particle per cell
        // sim.Lx = sim.dx;
        // sim.particle_volume = sim.dx*sim.dx;
        // sim.particle_mass = sim.rho * sim.particle_volume;
        // sim.particles = Particles(sim.Np);
        // for(int p = 0; p < sim.Np; p++){
        //     TV vec = TV::Ones()*0.5*sim.dx;
        //     vec(1) = (0.5+p)*sim.dx;
        //     sim.particles.x[p] = vec;
        // }
    }

    else if (setup == 1){     ///// Gran collapse:
        sim.pbc = false;
        sim.Lx = 0.08;
        sim.Ly = 0.05;
        #ifdef THREEDIM
            sim.Lz = 0.02;
            // SampleParticles(sim.Lx, sim.Ly, sim.Lz,  0.0004, sim);
            SampleParticles(sim.Lx, sim.Ly, sim.Lz,  0.0003, sim);
        #else
            // SampleParticles(sim.Lx, sim.Ly,    0.0003, 6, 0, sim);
            // SampleParticles(sim.Lx, sim.Ly,    0.0001, 6, 0, sim);
            SampleParticles(sim.Lx, sim.Ly,    0.00013, 6, 0, sim);
        #endif
        for(int p = 0; p < sim.Np; p++){
            sim.particles.x[p](1) += 0.5*sim.dx;
        }
        auto new_part_x = sim.particles.x;
        #ifdef THREEDIM
            TV tmp_part(sim.Lx/2.0, 0.0, sim.Lz/2.0);
        #else
            TV tmp_part(sim.Lx/2.0, 0.0);
        #endif
        new_part_x.push_back(tmp_part);
        sim.Np += 1;
        sim.particles = Particles(sim.Np);
        sim.particles.x = new_part_x;
    }

    else if (setup == 2){     ///// Incl slope:
        sim.pbc = false;
        sim.Lx = 0.20;
        sim.Ly = 0.14;
        #ifdef THREEDIM
            sim.Lz = 0.1;
            // SampleParticles(sim.Lx, sim.Ly, sim.Lz,  0.0013, sim);
            SampleParticles(sim.Lx, sim.Ly, sim.Lz,  0.001, sim);
        #else
            SampleParticles(sim.Lx, sim.Ly,     0.001, 6, 0, sim); // comp3d
            // SampleParticles(sim.Lx, sim.Ly,    0.0006, 6, 0, sim); // dx2
            // SampleParticles(sim.Lx, sim.Ly,    0.00035, 6, 0, sim); // dx3
        #endif
        for(int p = 0; p < sim.Np; p++){
            sim.particles.x[p](1) += 0.5*sim.dx;
        }
        auto new_part_x = sim.particles.x;
        #ifdef THREEDIM
            TV tmp_part(sim.Lx/2.0, 0.0, sim.Lz/2.0);
        #else
            TV tmp_part(sim.Lx/2.0, 0.0);
        #endif
        new_part_x.push_back(tmp_part);
        sim.Np += 1;
        sim.particles = Particles(sim.Np);
        sim.particles.x = new_part_x;
    }
    else {                    ///// Feeder with ramp
        sim.pbc = false;
        sim.Lx = 1.8;
        sim.Ly = 0.1;
        h_gate = 0.05;
        l_gate = 0.1;
        #ifndef THREEDIM
            // SampleParticles(sim.Lx, sim.Ly,     0.002, 6, 0, sim);
            // SampleParticles(sim.Lx, sim.Ly,     0.0005, 6, 0, sim); // dx3
            SampleParticles(sim.Lx, sim.Ly,     0.0004, 6, 0, sim); // dx35
            // SampleParticles(sim.Lx, sim.Ly,     0.0003, 6, 0, sim); // dx4
        #endif
        for(int p = 0; p < sim.Np; p++){
            sim.particles.x[p](0) -= sim.Lx;
            sim.particles.x[p](1) += l_gate + 0.5*sim.dx;
        }
        auto new_part_x = sim.particles.x;
        #ifndef THREEDIM
            TV tmp_part(l_gate, 0.0);
            new_part_x.push_back(tmp_part);
        #endif
        sim.Np += 1;
        sim.particles = Particles(sim.Np);
        sim.particles.x = new_part_x;
    }


    sim.dt_max = 0.5 * sim.dx / sim.wave_speed;
    // sim.dt_max = 0.5 * sim.dx / (std::sqrt(1.0e7/sim.rho));

    std::string name;

    if (setup == 0){    //// PBC
        #ifndef THREEDIM
            name = "Ground"; InfinitePlate ground = InfinitePlate(0,  1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,     1,0);  sim.objects.push_back(ground);
        #endif
    }
    else if (setup == 1){
        #ifdef THREEDIM //// 3D Grancular collapse
            name = "Ground";     InfinitePlate ground     = InfinitePlate(0,      1e20, -1e20, bottom, STICKY,   0.0, name,   0,0,   0,   1,0);  sim.objects.push_back(ground);
            name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      1e20,     0, left,   SEPARATE, 0.0, name,   0,0.45,0,   1,0);  sim.objects.push_back(side_left);
            name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx, 1e20,     0, right,  SEPARATE, 0.0, name,   0,0.45,0,   1,0);  sim.objects.push_back(side_right);
            name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0,      1e20, -1e20, back,   SEPARATE, 0.0, name,   0,0,   0,   1,0);  sim.objects.push_back(side_back);
            name = "SideFront";  InfinitePlate side_front = InfinitePlate(sim.Lz, 1e20, -1e20, front,  SEPARATE, 0.0, name,   0,0,   0,   1,0);  sim.objects.push_back(side_front);
        #else          //// 2D Granular collapse
            name = "Ground";     InfinitePlate ground     = InfinitePlate(0,       1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,        1,0);  sim.objects.push_back(ground);
            name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,       1e20, 0,     left,   SEPARATE, 0.0, name,   0, 0.45,     1,0);  sim.objects.push_back(side_left);
            name = "SideRight";  InfinitePlate side_right = InfinitePlate(sim.Lx,  1e20, 0,     right,  SEPARATE, 0.0, name,   0, 0.45,     1,0);  sim.objects.push_back(side_right);
        #endif
    }
    else if (setup == 2){
        #ifdef THREEDIM //// 3D Inclined slope
            name = "Ground";     InfinitePlate ground     = InfinitePlate(0,      1e20, -1e20, bottom,  STICKY,  0.0, name,  0, 0, 0,   1,0);  sim.objects.push_back(ground);
            name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      1e20, -1e20, left,   SEPARATE, 0.0, name,  0, 0, 0,   1,0);  sim.objects.push_back(side_left);
            name = "SideBack";   InfinitePlate side_back  = InfinitePlate(0,      1e20, -1e20, back,   SEPARATE, 0.0, name,  0, 0, 0,   1,0);  sim.objects.push_back(side_back);
            name = "SideFront";  InfinitePlate side_front = InfinitePlate(sim.Lz, 1e20, -1e20, front,  SEPARATE, 0.0, name,  0, 0, 0,   1,0);  sim.objects.push_back(side_front);
        #else           //// 2D Inclined slope
            name = "Ground";     InfinitePlate ground     = InfinitePlate(0,       1e20, -1e20, bottom,  STICKY,   0.0, name,  0, 0,   1,0);  sim.objects.push_back(ground);
            name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,       1e20, -1e20,  left,   SEPARATE, 0.0, name,  0, 0,   1,0);  sim.objects.push_back(side_left);
        #endif
    }
    else {
        #ifndef THREEDIM
            //// 2D Feeder
            name = "Ground";     InfinitePlate ground     = InfinitePlate(0,      1e20, -1e20, bottom, STICKY,   0.0, name,   0, 0,     1,0);  sim.objects.push_back(ground);
            name = "SideRight";  InfinitePlate side_right = InfinitePlate(l_gate, 1e20, h_gate, right, SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_right);
            // name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      1e20, -1e20, left,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_left);
            // name = "Top";        InfinitePlate lid        = InfinitePlate(h_gate, 1e20, l_gate, top,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(lid);
            name = "Ramp";       InfinitePlate ramp       = InfinitePlate(l_gate,      0, -1e20, bottom, SEPARATE, 0.2, name,   0, 0,     1,0);  sim.objects.push_back(ramp);
            name = "SideLeft";   InfinitePlate side_left  = InfinitePlate(0,      l_gate, -1e20, left,   SEPARATE, 0.0, name,   0, 0,     1,0);  sim.objects.push_back(side_left);
        #endif
    }

    // Elastoplasticity
    sim.elastic_model = StvkWithHencky;
    // sim.plastic_model = NoPlasticity;
    // sim.plastic_model = VonMises;
    // sim.plastic_model = DruckerPrager;
    // sim.plastic_model = PerzynaVM;
    // sim.plastic_model = PerzynaDP;
    // sim.plastic_model = PerzynaMuIDP;
    // sim.plastic_model = MCC;
    // sim.plastic_model = MCCHard;
    // sim.plastic_model = MCCHardExp;
    // sim.plastic_model = PerzynaMCC;
    // sim.plastic_model = PerzynaMCCHard;
    sim.plastic_model = PerzynaMuIMCC;
    // sim.plastic_model = PerzynaSinterMCC;

    sim.use_jop_definitions = true;

    if (setup == 2) // incl slope
        sim.dp_slope = std::tan(24.5 * M_PI / 180);
    else // feeder and gran coll
        sim.dp_slope = std::tan(20.9 * M_PI / 180);

    sim.dp_cohesion = 0;

    sim.vm_ptensile = -5e10;
    sim.yield_stress_orig = 5e3;
    sim.yield_stress_min = sim.yield_stress_orig;

    sim.perzyna_exp = 1;
    sim.perzyna_visc = 1;

    sim.beta = 0.0;
    sim.M = sim.dp_slope;
    sim.p0 = 1e2;

    sim.xi = 50; // In case of sintering, xi = 1/(lambda*phi_0)

    sim.xi_nonloc = 0;
    sim.nonlocal_l = 0;

    // For mu(I) rheology only
    if (setup == 0)
        sim.grain_diameter = 2.5e-3;
    else if (setup == 1)
        sim.grain_diameter = 2e-3; // gran collapse
    else if (setup == 2)
        sim.grain_diameter = 0.7e-3; // incl slope
    else
        sim.grain_diameter = 2.5e-3; // feeder

    sim.rho_s           = 2500;
    sim.in_numb_ref     = 0.279;
    sim.mu_1            = sim.dp_slope; //sim.M
    sim.mu_2            = std::tan(32.76 * M_PI / 180);

    //// For sintering only
    // sim.sinter_tc = 20;   // 1     // 20
    // sim.sinter_ec = 0.03; // 1e-3  // 0.03
    // sim.sinter_Sc = 20;
    // sim.particles.sinter_S.resize(sim.Np); std::fill( sim.particles.sinter_S.begin(), sim.particles.sinter_S.end(), sim.sinter_Sinf );

    // For MCC only:
    T eps_pl_vol_init = -std::asinh(sim.p0/sim.K) / sim.xi;
    sim.particles.eps_pl_vol_mcc.resize(sim.Np); std::fill( sim.particles.eps_pl_vol_mcc.begin(), sim.particles.eps_pl_vol_mcc.end(), eps_pl_vol_init );

    sim.simulate();
    // sim.validateRMA();


	return 0;
}
