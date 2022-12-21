#ifndef SAMPLING_PARTICLES_HPP
#define SAMPLING_PARTICLES_HPP

#include "tools.hpp"
#include "data_structures.hpp"
#include "poisson_disk_sampling.hpp"

#ifdef THREEDIM

    template <typename S>
    void SampleParticles(const T Lx, const T Ly, const T Lz, T kRadius, S& sim){
        std::uint32_t kAttempts = 30;
        std::uint32_t kSeed = 42;
        std::array<T, 3> kXMin = std::array<T, 3>{{0, 0, 0}};
        std::array<T, 3> kXMax = std::array<T, 3>{{Lx, Ly, Lz}};

        debug("Sampling particles...");
        std::vector<std::array<T, 3>> samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);
        for(int p = 0; p < sim.Np; p++){
            for(int d = 0; d < 3; d++){
                sim.particles.x[p](d) = samples[p][d];
            }
        }

        unsigned int Npx = Lx / (Lx*Ly*Lz) * std::pow( std::pow(Lx*Ly*Lz, 2) * sim.Np, 1.0/3.0);
        T dx_p = (Lx / Npx);
        sim.dx = 2 * dx_p;
        sim.particle_volume = std::pow(dx_p, 3);
        sim.particle_mass = sim.rho * sim.particle_volume;
    } // end SampleParticles

#else // TWODIM

    template <typename S>
    void SampleParticles(const T Lx, const T Ly, T kRadius, T ppc, unsigned int front_type, S& sim){
        std::uint32_t kAttempts = 30;
        std::uint32_t kSeed = 42;
        std::array<T, 2> kXMin = std::array<T, 2>{{0, 0}};
        std::array<T, 2> kXMax = std::array<T, 2>{{Lx, Ly}};

        debug("Sampling particles...");
        std::vector<std::array<T, 2>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        std::vector<std::array<T, 2>> samples;

        debug("    Number of square samples: ", square_samples.size());

        /////// Triangle
        if (front_type == 1){
            for(int p = 0; p < square_samples.size(); p++){
                if (square_samples[p][1] < Ly - Ly/Lx*square_samples[p][0]){
                    samples.push_back(square_samples[p]);
                }
            }
        }
        /////// Rounded Edge
        else if (front_type == 2){
            for(int p = 0; p < square_samples.size(); p++){
                if (square_samples[p][0] < Lx-Ly || square_samples[p][1] < std::sqrt(Ly*Ly - (square_samples[p][0]-Lx+Ly)*(square_samples[p][0]-Lx+Ly))){
                    samples.push_back(square_samples[p]);
                }
            }
        }

        /////// Quadratic Edge
        else if (front_type == 3){
            T ym = Ly;
            T xm = Lx - std::sqrt(ym);
            for(int p = 0; p < square_samples.size(); p++){
                if (square_samples[p][0] < xm || square_samples[p][1] < ym - std::pow(square_samples[p][0] - xm, 2) ){
                    samples.push_back(square_samples[p]);
                }
            }
        }

        /////// Double Quadratic Edge
        else if (front_type == 4){
            T ym  = Ly;
            T xmr = Lx - std::sqrt(ym);
            T xml = std::sqrt(ym);
            T xmm = Lx / 2.0;

            for(int p = 0; p < square_samples.size(); p++){

                if (square_samples[p][0] > xmm){
                    if (square_samples[p][0] < xmr || square_samples[p][1] < ym - std::pow(square_samples[p][0] - xmr, 2) ){
                        samples.push_back(square_samples[p]);
                    }
                } else{
                    if (square_samples[p][0] > xml || square_samples[p][1] < ym - std::pow(xml - square_samples[p][0], 2) ){
                        samples.push_back(square_samples[p]);
                    }
                }

            }
        }
        else{
            debug("    No front type specified (1,2,3), using just a square.");
            samples = square_samples;
        }

        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);
        for(int p = 0; p < sim.Np; p++){
            for(int d = 0; d < 2; d++){
                sim.particles.x[p](d) = samples[p][d];
            }
        }

        sim.dx = std::sqrt(ppc / T(square_samples.size()) * Lx*Ly);
        sim.particle_volume = sim.dx * sim.dx / ppc;
        sim.particle_mass = sim.rho * sim.particle_volume;

        debug("    dx set to ", sim.dx);
    } // end SampleParticles

#endif // DIMENSION


#endif  // SAMPLING_PARTICLES_HPP
