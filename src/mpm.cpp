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
//  * make Nx and Ny for grid

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
inline double wip(TV2 xp, double xi, double yi, double h){
    return N( (xp(0) - xi) / h ) * N( (xp(1) - yi) / h );
}

TV2 grad_wip(TV2 xp, double xi, double yi, double h){
    TV2 out;
    out(0) = 0;
    out(1) = 0;
    return out;
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
  TV2 F;
};


class Simulation{
public:
	std::vector<Particle> particles;
  unsigned int Np;
  double particle_mass;
  double particle_volume;
  double rho;
  double E;
  double nu;
	double dt;
	double dx;

  double mu;
  double lambda;
  void getElasticParams(){
      mu = 0;
      lambda = 0;
  }

  unsigned int N;
  Eigen::MatrixXd grid_X;
  Eigen::MatrixXd grid_Y;
  Eigen::MatrixXd grid_VX;
  Eigen::MatrixXd grid_VY;
  Eigen::MatrixXd grid_mass;



	void advanceStep();

	// advanceStep uses the following:
	void P2G();
  void explicitEulerUpdate();
	void G2P();

};

void Simulation::advanceStep(){

}

void Simulation::P2G(){

}

void Simulation::G2P(){

}

void Simulation::explicitEulerUpdate(){
    TV2 force, velocity_increment;
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            force = Eigen::Vector2d::Zero();
            for(int p=0; p<Np; p++){
                Eigen::JacobiSVD<Eigen::Matrix2d> svd((particles[p].F * particles[p].F.transpose()).array().log(), Eigen::ComputeThinU | Eigen::ComputeThinV);
                TV2 sigma = svd.singularValues().array().exp().sqrt();
                TV2 logSigma = sigma.array().log();
                TV2 invSigma = sigma.inverse();
                TV2 dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.sum()*invSigma ) * svd.matrixV().transpose();
                force -= particle_volume * dPsidF * particles[p].F.transpose() * grad_wip(particles[p].x, grid_X(i,j), grid_Y(i,j), dx);
            }
            velocity_increment = dt * force / grid_mass(i,j);
            grid_VX(i,j) += velocity_increment(0);
            grid_VY(i,j) += velocity_increment(1);
        }
    }
}


int main(){
    std::cout << "This is Larsie" << std::endl;
		return 0;
}
