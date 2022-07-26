#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include "poisson_disk_sampling.hpp"


//// PARAMETERS ////
// typedef float T;
typedef double T;
#define CUBICSPLINES
// #define THREEDIM // Uncomment for 2D
// #define DIMENSION 3 // Needed for OMP collapse
#define DIMENSION 2 // Needed for OMP collapse
#define TINYPLY_IMPLEMENTATION // Use tinyply

// #define WARNINGS // if write warnings to screen
////////////////////////

#ifdef THREEDIM
    typedef Eigen::Matrix<T, 3, 3> TM;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
    typedef Eigen::Matrix<T, 3, 1> TV;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
    typedef Eigen::Array<T,3,1> TA;
#else
    typedef Eigen::Matrix<T, 2, 2> TM;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
    typedef Eigen::Matrix<T, 2, 1> TV;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
    typedef Eigen::Array<T,2,1> TA;
#endif
////////////////////////

enum PlateType { top, bottom, left, right, front, back };
enum ElasticModel { StvkWithHencky, NeoHookean };
enum PlasticModel { NoPlasticity, VonMises, DruckerPrager, DPSoft, ModifiedCamClay, ModifiedCamClayHard, PerzynaVM, PerzynaDP, PerzynaMuIDP, PerzynaMCC, PerzynaMuIMCC, PerzynaSinterMCC};
enum BoundaryCondition { STICKY, SLIP, SEPARATE };

///////////////////// TOOLS ////////////////////////

template <typename T>
void debug(T in){
  std::cout << in << std::endl;
}
template <typename T, typename U>
void debug(T in1, U in2){
  std::cout << in1 << in2 << std::endl;
}
template <typename T, typename U, typename V>
void debug(T in1, U in2, V in3){
  std::cout << in1 << in2 << in3 << std::endl;
}
template <typename T, typename U, typename V, typename W>
void debug(T in1, U in2, V in3, W in4){
  std::cout << in1 << in2 << in3 << in4 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X>
void debug(T in1, U in2, V in3, W in4, X in5){
  std::cout << in1 << in2 << in3 << in4 << in5 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y, typename Z>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6, Z in7){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << in7 << std::endl;
}

inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline T selfDoubleDot(TM& A){

    #ifdef THREEDIM
    T out = A(0,0)*A(0,0) + A(0,1)*A(0,1) + A(0,2)*A(0,2)
          + A(1,0)*A(1,0) + A(1,1)*A(1,1) + A(1,2)*A(1,2)
          + A(2,0)*A(2,0) + A(2,1)*A(2,1) + A(2,2)*A(2,2);
    #else
    T out = A(0,0)*A(0,0) + A(0,1)*A(0,1)
          + A(1,0)*A(1,0) + A(1,1)*A(1,1);
    #endif

    return out;
}

unsigned int load_array(std::vector<TV>& array, std::string file_name);
std::vector<T> linspace(T a, T b, size_t N);
bool copy_file(std::string source, std::string destination);

bool ModifiedCamClayHardRMA(T& p, T& q, int& exit, T M, T epv, T beta, T mu, T K, T xi);
bool ModifiedCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K);
bool PerzynaMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc);
bool PerzynaSinterMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T epv, T S, T visc, T Sinf, T tc, T ec, T plasitic_index);

// DO NOT USE THESE:
bool PerzynaCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc);
bool PerzynaQuadRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc);
bool CamClayRMA(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus);
bool QuadRMA(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus);

#ifdef CUBICSPLINES

    inline T N(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return 0.5 * uabs*uabs*uabs - uabs*uabs + 0.6666666666666666666666666666666666666666666666666;
        }
            else if (uabs < 2.0){
            return 0.1666666666666666666666666666666666666666666 * (2.0 - uabs) * (2.0 - uabs) * (2.0 - uabs);
        }
            else {
            return 0;
        }
    }

    inline T dNdu(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return u * (1.5 * uabs - 2.0);
        }
        else if (uabs < 2.0){
            return -0.5 * sgn(u) * (2.0 - uabs) * (2.0 - uabs);
        }
        else {
            return 0;
        }
    }

    inline T d2Ndu2(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return (3.0 * uabs - 2.0);
        }
        else if (uabs < 2.0){
            return (2.0 - uabs);
        }
        else {
            return 0;
        }
    }

#else // QUADRATIC SPLINES

    inline T N(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return 0.75 - uabs * uabs;
        }
    		else if (uabs < 1.5){
            return 0.5 * (1.5 - uabs) * (1.5 - uabs);
        }
    		else {
            return 0;
    	}
    }

    inline T dNdu(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return (-2*u);
        }
    	else if (uabs < 1.5){
            return (u - 1.5*sgn(u));
        }
    	else {
            return 0;
    	}
    }

    inline T d2Ndu2(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return -2.0;
        }
    	else if (uabs < 1.5){
            return 1;
        }
    	else {
            return 0;
    	}
    }

#endif


#ifdef THREEDIM

    inline T wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h){
        return N( (xp - xi) * one_over_h ) * N( (yp - yi) * one_over_h ) * N( (zp - zi) * one_over_h );
    }

    inline TV grad_wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h){
        TV out;
        out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * N((zp - zi) * one_over_h) * one_over_h,
               dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * N((zp - zi) * one_over_h) * one_over_h,
               dNdu((zp - zi) * one_over_h) * N((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h;
        return out;
    }

    inline T laplace_wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h, T one_over_h_square){
        T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) *  N((zp - zi) * one_over_h);
        T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) *  N((zp - zi) * one_over_h);
        T term3 = d2Ndu2((zp - zi) * one_over_h) * N((xp - xi) * one_over_h) *  N((yp - yi) * one_over_h);
        return ( term1 + term2 + term3 ) * one_over_h_square;
    }

#else // TWODIM

    inline T wip(T xp, T yp, T xi, T yi, T one_over_h){
        return N( (xp - xi) * one_over_h ) * N( (yp - yi) * one_over_h );
    }

    inline TV grad_wip(T xp, T yp, T xi, T yi, T one_over_h){
        TV out;
        out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h,
               dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * one_over_h;
        return out;
    }

    // inline T wip(T xp, T yp, T xi, T yi, T one_over_h){
    //     T Lx = 0.1;
    //     return N( fmod(xp - xi, Lx) * one_over_h ) * N( (yp - yi) * one_over_h );
    // }
    // inline TV grad_wip(T xp, T yp, T xi, T yi, T one_over_h){
    //     T Lx = 0.1;
    //     T xpminxi = xp-xi;
    //     if (xpminxi <= -Lx || xpminxi >= Lx)
    //         xpminxi = -fmod(xpminxi, Lx);
    //     TV out;
    //     out << dNdu( xpminxi  * one_over_h) * N((yp - yi) * one_over_h) * one_over_h,
    //            dNdu((yp - yi) * one_over_h) * N( xpminxi  * one_over_h) * one_over_h;
    //     return out;
    // }

    inline T laplace_wip(T xp, T yp, T xi, T yi, T one_over_h, T one_over_h_square){
        T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h);
        T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h);
        return ( term1 + term2 ) * one_over_h_square;
    }

#endif



class Particles{
public:
  Particles(unsigned int Np = 1){
      x.resize(Np);    std::fill( x.begin(),    x.end(),    TV::Zero() );
      v.resize(Np);    std::fill( v.begin(),    v.end(),    TV::Zero() );
      pic.resize(Np);  std::fill( pic.begin(),  pic.end(),  TV::Zero() );
      flip.resize(Np); std::fill( flip.begin(), flip.end(), TV::Zero() );

      eps_pl_dev.resize(Np);     std::fill( eps_pl_dev.begin(),     eps_pl_dev.end(),     0.0 );
      eps_pl_vol.resize(Np);     std::fill( eps_pl_vol.begin(),     eps_pl_vol.end(),     0.0 );
      eps_pl_vol_mcc.resize(Np);     std::fill( eps_pl_vol_mcc.begin(),     eps_pl_vol_mcc.end(),     0.0 );
      eps_pl_vol_pradhana.resize(Np);std::fill( eps_pl_vol_pradhana.begin(),eps_pl_vol_pradhana.end(),0.0 );

      // eps_pl_vol_abs.resize(Np);     std::fill( eps_pl_vol_abs.begin(),     eps_pl_vol_abs.end(),     0.0 );
      // fail_crit.resize(Np); std::fill( fail_crit.begin(), fail_crit.end(), false );

      eps_pl_dev_nonloc.resize(Np);  std::fill( eps_pl_dev_nonloc.begin(),  eps_pl_dev_nonloc.end(),  0.0 );
      delta_gamma_nonloc.resize(Np); std::fill( delta_gamma_nonloc.begin(), delta_gamma_nonloc.end(), 0.0 );
      delta_gamma.resize(Np);        std::fill( delta_gamma.begin(),        delta_gamma.end(),        0.0 );

      sinter_S.resize(Np); std::fill( sinter_S.begin(), sinter_S.end(), 0.0 );

      yield_stress_orig.resize(Np); std::fill( yield_stress_orig.begin(), yield_stress_orig.end(), 0.0 );

      viscosity.resize(Np); std::fill( viscosity.begin(), viscosity.end(), 0.0 );
      muI.resize(Np); std::fill( muI.begin(), muI.end(), 0.0 );

      hencky.resize(Np);             std::fill( hencky.begin(),             hencky.end(),      TV::Zero() );

      tau.resize(Np); std::fill( tau.begin(), tau.end(), TM::Zero()     );
      F.resize(Np);   std::fill( F.begin(),   F.end(),   TM::Identity() );
  }

  std::vector<TV> x;
  std::vector<TV> x0;
  std::vector<TV> v;
  std::vector<TV> pic;
  std::vector<TV> flip;

  std::vector<T> eps_pl_dev;
  std::vector<T> eps_pl_vol;
  std::vector<T> eps_pl_vol_mcc;
  std::vector<T> eps_pl_vol_pradhana;
  // std::vector<T> eps_pl_vol_abs;
  // std::vector<bool> fail_crit;
  std::vector<T> eps_pl_dev_nonloc;
  std::vector<T> delta_gamma_nonloc;
  std::vector<T> delta_gamma;

  std::vector<T> sinter_S;

  std::vector<T> yield_stress_orig;

  std::vector<T> viscosity;
  std::vector<T> muI;

  std::vector<TV> hencky;

  std::vector<T> cohesion_proj;

  std::vector<TM> tau;
  std::vector<TM> F;

};

class Grid{
public:
    Grid(){}
    std::vector<T> x;
    std::vector<T> y;
    #ifdef THREEDIM
    std::vector<T> z;
    #endif
    std::vector<TV> v;
    std::vector<TV> flip;
    std::vector<T> mass;
    std::vector<T> delta_gamma;
    T xc;
    T yc;
    #ifdef THREEDIM
    T zc;
    #endif
};



#ifdef THREEDIM
    template <typename S>
    void SampleInBox(const T Lx, const T Ly, const T Lz, T kRadius, S& sim){
        std::uint32_t kAttempts = 30;
        std::uint32_t kSeed = 42;
        std::array<T, 3> kXMin = std::array<T, 3>{{0, 0, 0}};
        std::array<T, 3> kXMax = std::array<T, 3>{{Lx, Ly, Lz}};

        debug("Sampling particles...");
        std::vector<std::array<T, 3>> samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        sim.Np = samples.size();
        debug("Number of particles samples: ", sim.Np);

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
    } // end SampleInBox

#else // TWODIM

template <typename S>
void SampleInBox(const T Lx, const T Ly, T kRadius, S& sim){
    std::uint32_t kAttempts = 30;
    std::uint32_t kSeed = 42;
    std::array<T, 2> kXMin = std::array<T, 2>{{0, 0}};
    std::array<T, 2> kXMax = std::array<T, 2>{{Lx, Ly}};

    debug("Sampling particles...");
    std::vector<std::array<T, 2>> samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
    sim.Np = samples.size();
    debug("Number of particles samples: ", sim.Np);

    sim.particles = Particles(sim.Np);
    for(int p = 0; p < sim.Np; p++){
        for(int d = 0; d < 2; d++){
            sim.particles.x[p](d) = samples[p][d];
        }
    }

    unsigned int Npx = std::sqrt(Lx/Ly * sim.Np);
    T dx_p = (Lx / Npx);
    sim.dx = 2 * dx_p;
    sim.particle_volume = std::pow(dx_p, 2);
    sim.particle_mass = sim.rho * sim.particle_volume;
} // end SampleInBox

#endif // DIMENSION

template <typename S>
void SampleIn2DSpecial(const T Lx, const T Ly, T kRadius, double ppc, unsigned int front_type, S& sim){
    std::uint32_t kAttempts = 30;
    std::uint32_t kSeed = 42;
    std::array<T, 2> kXMin = std::array<T, 2>{{0, 0}};
    std::array<T, 2> kXMax = std::array<T, 2>{{Lx, Ly}};

    debug("Sampling particles...");
    std::vector<std::array<T, 2>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
    std::vector<std::array<T, 2>> samples;

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
        for(int p = 0; p < square_samples.size(); p++){
            T ym = Ly;
            T xm = Lx - std::sqrt(ym);
            if (square_samples[p][0] < xm || square_samples[p][1] < ym - std::pow(square_samples[p][0] - xm, 2) ){
                samples.push_back(square_samples[p]);
            }
        }
    }
    else{
        debug("No front type specified (1,2,3), using just a square.");
        samples = square_samples;
    }

    sim.Np = samples.size();
    debug("Number of particles samples: ", sim.Np);

    sim.particles = Particles(sim.Np);
    for(int p = 0; p < sim.Np; p++){
        for(int d = 0; d < 2; d++){
            sim.particles.x[p](d) = samples[p][d];
        }
    }


    unsigned int Npx = std::sqrt(Lx/Ly * sim.Np);
    T dx_p = (Lx / Npx);
    sim.dx = std::sqrt(ppc) * dx_p;
    sim.particle_volume = std::pow(dx_p, 2);
    sim.particle_mass = sim.rho * sim.particle_volume;
} // end SampleIn2DTriangle




#endif  // TOOLS_HPP
