#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

// TODO:
// 	* Fix type template for double/float
//  * 3D
//  * Cubic spline
//  * adaptive timestepping
//  * complete doxygen

typedef Eigen::Vector2d TV2;
typedef Eigen::Matrix2d TM2;

///////////////////// TOOLS ////////////////////////

void debug(std::string in){
  std::cout << in << std::endl;
}
void debug(int in){
  std::cout << in << std::endl;
}
void debug(double in){
  std::cout << in << std::endl;
}
void debug(std::string in1, int in2){
  std::cout << in1 << in2 << std::endl;
}

/*!
 \param x x
 \return weight

 Quadratic spline basis function
*/
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
/*!
 \param u u
 \return derivative of function evaluated at u

 Derivative of quadratic spline basis function
*/
inline double dNdu(double u){
    double uabs = std::abs(u);
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
inline double wip(double xp, double yp, double xi, double yi, double h){
    return N( (xp - xi) / h ) * N( (yp - yi) / h );
}

inline double gradx_wip(double xp, double yp, double xi, double yi, double h){
    return ( dNdu((xp - xi) / h) *  N((yp - yi) / h) ) / h;
}
inline double grady_wip(double xp, double yp, double xi, double yi, double h){
    return ( dNdu((yp - yi) / h) *  N((xp - xi) / h) ) / h;
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
  Particle(){
      F = TM2::Identity();
  }
  TM2 F;
};

class Simulation{
public:
  Simulation();

  void setElasticParams(double E, double nu, double density){
      mu = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) );
      lambda = E / (2.0*(1.0+nu));
      rho = density;
      dt = 0.1 * dx * std::sqrt(rho/E);
  }

  void advanceStep(){
      remesh();
      debug("Finished remesh");
      P2G();
      debug("Finished P2G");
      explicitEulerUpdate();
      debug("Finished explicitEulerUpdate");
      G2P();
      debug("Finished G2P");
  };

  void simulate(){
      double t = 0;
      for (int i = 0; i < Nt; i++){
          advanceStep();
          t += dt;
          std::cout << "Step: " << i << "\t Time: " << t << std::endl;
      }
      saveSim();
  };

  void saveSim(){
    std::ofstream outFile("out.csv");
    for(int p = 0; p < Np; p++){
        outFile << p << "," << particles_x[p] << "\t" << particles_y[p] << "\t" << particles_vx[p] << "\t" << particles_vy[p] << "\n";
    }
  };

private:

  unsigned int Nt;
  double dt;
  double dx;

  double rho;
  double mu;
  double lambda;

  unsigned int Np;
  double particle_mass;
  double particle_volume; // initial particle volume V0

	std::vector<Particle> particles;
  Eigen::VectorXd particles_x;
  Eigen::VectorXd particles_y;
  Eigen::VectorXd particles_vx;
  Eigen::VectorXd particles_vy;

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
    for(int i=0; i<Nx; i++){
        //debug("i = ", i);
        for(int j=0; j<Ny; j++){
            //debug("j = ", j);
            double xi = grid_X(i,j);
            double yi = grid_Y(i,j);
            double vxi = 0;
            double vyi = 0;
            double mass = 0;
            for(int p=0; p<Np; p++){
                //debug(p);
                double xp = particles_x(p);
                double yp = particles_y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
                    double weight = wip(xp, yp, xi, yi, dx);
                    mass += weight;
                    vxi  += particles_vx(p) * weight;
                    vyi  += particles_vy(p) * weight;
                }
            }
            //debug(vxi);
            grid_mass(i,j) = mass * particle_mass;
            grid_VX(i,j)   = vxi  * particle_mass;
            grid_VY(i,j)   = vyi  * particle_mass;
        }
    }
}

void Simulation::G2P(){
    for(int p=0; p<Np; p++){
        double xp = particles_x(p);
        double yp = particles_y(p);
        double vxp = 0;
        double vyp = 0;
        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){
                double xi = grid_X(i,j);
                double yi = grid_Y(i,j);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
                    double weight = wip(xp, yp, xi, yi, dx);
                    vxp += grid_VX(i,j) * weight;
                    vyp += grid_VY(i,j) * weight;
                }
            }
        }
        particles_vx(p) = vxp;
        particles_vy(p) = vyp;
    }
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

    //std::cout << grid_X.rows() << std::endl;
    //std::cout << grid_X.cols() << std::endl;
    grid_VX   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_VY   = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());
    grid_mass = Eigen::MatrixXd::Zero(grid_X.rows(), grid_X.cols());

}

void Simulation::explicitEulerUpdate(){
    TV2 grad_wip, sigma;
    TM2 Fe, logSigma, invSigma, dPsidF;
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            double xi = grid_X(i,j);
            double yi = grid_Y(i,j);
            TV2 force = TV2::Zero();
            for(int p=0; p<Np; p++){
                double xp = particles_x(p);
                double yp = particles_y(p);
                if ( std::abs(xp-xi) < 1.5*dx && std::abs(xp-yi) < 1.5*dx){
                    Fe = particles[p].F; // must update F!!
                    Eigen::JacobiSVD<TM2> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
                    logSigma = sigma.array().abs().max(1e-4).log().matrix().asDiagonal();
                    invSigma = sigma.array().abs().max(1e-4).matrix().asDiagonal().inverse();
                    dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
                    std::cout << invSigma << std::endl << std::endl;
                    grad_wip(0) = gradx_wip(xp, yp, xi, yi, dx);
                    grad_wip(1) = grady_wip(xp, yp, xi, yi, dx);
                    force -= particle_volume * dPsidF * Fe.transpose() * grad_wip;
                }
            }
            TV2 velocity_increment = dt * force / grid_mass(i,j);
            grid_VX(i,j) += velocity_increment(0);
            grid_VY(i,j) += velocity_increment(1);
        }
    }
}

Simulation::Simulation(){
  // default case: unit box with 10 times dx = 0.1
    Nt = 1;
    dx = 0.1;
    rho = 250;
    setElasticParams(1e6, 0.3, rho);

    Nx = 10;
    Ny = 10;
    Np = Nx*Ny*4;
    debug("Np = ", Np);

    particle_volume = 1.0 / Np;
    particle_mass = rho * 1.0 / Np;

    grid_X = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_Y = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_VX = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_VY = Eigen::MatrixXd::Zero(Nx, Ny);
    grid_mass = Eigen::MatrixXd::Zero(Nx, Ny);

    Particle part;
    particles.resize(Np);
    std::fill(particles.begin(), particles.end(), part);

    particles_x  = Eigen::VectorXd::Zero(Np);
    particles_y  = Eigen::VectorXd::Zero(Np);
    particles_vx = Eigen::VectorXd::Zero(Np);
    particles_vy = Eigen::VectorXd::Zero(Np);

    int p = -1;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){

            p++;
            double px = (i+0.25)*dx;
            double py = (j+0.25)*dx;
            double pvx = 0.01*std::sin(M_PI*px);
            double pvy = 0.01*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.75)*dx;
            py = (j+0.75)*dx;
            pvx = 0.01*std::sin(M_PI*px);
            pvy = 0.01*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.25)*dx;
            py = (j+0.75)*dx;
            pvx = 0.01*std::sin(M_PI*px);
            pvy = 0.01*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;

            p++;
            px = (i+0.75)*dx;
            py = (j+0.25)*dx;
            pvx = 0.01*std::sin(M_PI*px);
            pvy = 0.01*std::sin(M_PI*py);
            particles_x(p) = px;
            particles_y(p) = py;
            particles_vx(p) = pvx;
            particles_vy(p) = pvy;
        } // end for i
    } // end for j
    debug("Final p = ", p);
    debug("Finished constructor");
}


int main(){
    std::cout << "This is Larsie" << std::endl;
    Simulation sim;
    sim.simulate();
		return 0;
}
