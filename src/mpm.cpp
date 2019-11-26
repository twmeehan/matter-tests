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
typedef Eigen::Matrix2d TM2;
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

TV2 grad_wip(double xp, double yp, double xi, double yi, double h){
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

class Particle{
public:
  TM2 F;
};

class Simulation{
public:
  double dt;
  double dx;

  unsigned int Np;
  double particle_mass;
  double particle_volume;

	std::vector<Particle> particles;
  Eigen::VectorXd particles_x;
  Eigen::VectorXd particles_y;
  Eigen::VectorXd particles_vx;
  Eigen::VectorXd particles_vy;

  double rho;
  double mu;
  double lambda;

  void setElasticParams(double E, double nu){
      mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
      lambda = E / (2.0*(1.0+nu));
  }

  void advanceStep(){
      remesh();
      P2G();
      explicitEulerUpdate();
      G2P();
  };

private:

  unsigned int Nx, Ny;
  Eigen::MatrixXd grid_X;
  Eigen::MatrixXd grid_Y;
  Eigen::MatrixXd grid_VX;
  Eigen::MatrixXd grid_VY;
  Eigen::MatrixXd grid_mass;

  void remesh();
	void P2G();
  void explicitEulerUpdate();
	void G2P();

}; // END mpm class

void Simulation::P2G(){

}

void Simulation::G2P(){

}

void Simulation::remesh(){
    double min_x = particles_x.minCoeff();
    double min_y = particles_y.minCoeff();
    double max_x = particles_x.maxCoeff();
    double max_y = particles_y.maxCoeff();

    double Lx = max_x - min_x;
    double Ly = max_y - min_y;

    Nx = std::ceil(Lx / dx);
    Ny = std::ceil(Ly / dx);

    double mid_x = min_x + Lx / 2.0;
    double mid_y = min_y + Ly / 2.0;

    double low_x  = mid_x - Nx*dx/2.0;
    double low_y  = mid_y - Ny*dx/2.0;
    double high_x = mid_x + Nx*dx/2.0;
    double high_y = mid_y + Ny*dx/2.0;

    Nx++;
    Ny++;

    Eigen::VectorXd lin_x = Eigen::VectorXd::LinSpaced(Nx, low_x, high_x);
    Eigen::VectorXd lin_y = Eigen::VectorXd::LinSpaced(Ny, low_y, high_y);

    meshgrid(lin_x, lin_y, grid_X, grid_Y);
}

void Simulation::explicitEulerUpdate(){
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            double xi = grid_X(i,j);
            double yi = grid_Y(i,j);
            TV2 force = TV2::Zero();
            for(int p=0; p<Np; p++){
                double xp = particles_x(p);
                double yp = particles_y(p);
                if ( std::abs(xp-xi) < 1.5*dx || std::abs(xp-yi) < 1.5*dx){
                    TM2 Fe = particles[p].F;
                    Eigen::JacobiSVD<TM2> svd((Fe*Fe.transpose()).array().log(), Eigen::ComputeThinU | Eigen::ComputeThinV);
                    TV2 sigma = svd.singularValues().array().exp().sqrt();
                    TM2 logSigma = sigma.array().log().matrix().asDiagonal();
                    TM2 invSigma = sigma.asDiagonal().inverse();
                    TM2 dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
                    force -= particle_volume * dPsidF * Fe.transpose() * grad_wip(xp, yp, xi, yi, dx);
                }
            }
            TV2 velocity_increment = dt * force / grid_mass(i,j);
            grid_VX(i,j) += velocity_increment(0);
            grid_VY(i,j) += velocity_increment(1);
        }
    }
}


int main(){
    std::cout << "This is Larsie" << std::endl;
		return 0;
}
