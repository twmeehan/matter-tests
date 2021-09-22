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

//// PARAMETERS ////
typedef float T;
#define CUBICSPLINES
// #define THREEDIM // Uncomment for 2D
// #define DIMENSION 3 // Needed for OMP collapse
#define DIMENSION 2 // Needed for OMP collapse

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

enum PlateType { top, bottom, left, right, front, back};
enum ElasticModel { StvkWithHencky, NeoHookean };
enum PlasticModel { NoPlasticity, VonMises, DPSimpleSoft, Curved, PerzynaVM, PerzynaDP };
enum BoundaryCondition { STICKY, SLIP };

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

bool PerzynaCamClayReturnMapping(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T perzyna_exp, T perzyna_visc);
bool CamClayReturnMapping(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus);
bool QuadraticReturnMapping(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus);
bool AnalQuadReturnMapping(T& p, T& q, int& exit, T M, T p0, T beta);

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

    inline T laplace_wip(T xp, T yp, T xi, T yi, T one_over_h, T one_over_h_square){
        T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h);
        T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h);
        return ( term1 + term2 ) * one_over_h_square;
    }

#endif

#endif  // TOOLS_HPP
