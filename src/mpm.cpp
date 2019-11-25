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

// quadratic spline
double N(double x){
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
double wip(TV2 xp, TV2 xi, double h){
    return N( (xp(0) - xi(0)) / h ) * N( (xp(1) - xi(1)) / h );
}

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


		return 0;
}
