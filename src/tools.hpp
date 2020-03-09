#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

/// FLOAT OR DOUBLE ////
typedef float T;
////////////////////////
typedef Eigen::Matrix<T, 2, 2> TM2;
typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
typedef Eigen::Matrix<T, 2, 1> TV2;
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
typedef Eigen::Array<T,2,1> TA2;


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

/*!
 \param x x
 \return weight

 Quadratic spline basis function
*/
inline T N(T x){
    T xabs = std::abs(x);
    if (xabs < 0.5){
        return 0.75 - xabs * xabs;
    }
		else if (xabs < 1.5){
        return 0.5 * (1.5 - xabs) * (1.5 - xabs);
    }
		else {
        return 0;
	}
}
/*!
 \param u u
 \return derivative of function evaluated at u

 Derivative of quadratic spline basis function
*/
inline T dNdu(T u){
    T uabs = std::abs(u);
    if (uabs < 0.5){
        return (-2*u);
    }
	else if (uabs < 1.5){
        return (-u);
    }
	else {
        return 0;
	}
}
inline T wip(T xp, T yp, T xi, T yi, T h){
    return N( (xp - xi) / h ) * N( (yp - yi) / h );
}

inline T gradx_wip(T xp, T yp, T xi, T yi, T h){
    return ( dNdu((xp - xi) / h) *  N((yp - yi) / h) ) / h;
}
inline T grady_wip(T xp, T yp, T xi, T yi, T h){
    return ( dNdu((yp - yi) / h) *  N((xp - xi) / h) ) / h;
}

#endif  // TOOLS_HPP
