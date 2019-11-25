#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

// TODO:
// 	* Fix type template for double/float
//  * 3D
//  * Cubic spline
//  * adaptive timestepping

typedef Eigen::Vector2d TV2;

///////////////////// TOOLS ////////////////////////

// quadratic spline
inline double N(double x){
    double xabs = std::abs(x);
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
inline double wip(TV2 xp, TV2 xi, double h){
    return N( (xp(0) - xi(0)) / h ) * N( (xp(1) - xi(1)) / h );
}

//  in : x, y column vectors, can be of type
//       X, Y matrices, used to save the mesh
void meshgrid(const Eigen::Matrix<double, -1, 1>& x,
              const Eigen::Matrix<double, -1, 1>& y,
              Eigen::Matrix<double, -1, -1>& X,
              Eigen::Matrix<double, -1, -1>& Y) {
  const long nx = x.size(), ny = y.size();
  X.resize(ny, nx);
  Y.resize(ny, nx);
  for (long i = 0; i < ny; ++i) {
    X.row(i) = x.transpose();
  }
  for (long j = 0; j < nx; ++j) {
    Y.col(j) = y;
  }
}

//////////////////////////////////////////////////////////////////

// move to particle.hpp
class Particle{
public:
	TV2 x;
	TV2 v;
	double mass;
};

class Simulation{
public:
	std::vector<Particle> particles;
	double dt;
	double dx;

	void advanceStep();

	// advanceStep uses the following:
	void P2G();
	void G2P();
	void explicitEuler();


};

void Simulation::advanceStep(){

}


int main(){
    std::cout << "This is Larsie" << std::endl;

		return 0;
}
