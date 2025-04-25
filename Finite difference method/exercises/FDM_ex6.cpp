#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>  // Required for Eigen compatibility
#include <cmath>
#include <iostream>
#include <fstream>

// PDE: ut = 0.1 uxx
// initial condition: u(x,0) = sin(pi*x)
// boundary conditions: u(0,t) = 0, u(1,t) = 0

namespace odeint = boost::numeric::odeint;

class HeatEquationSystem 
{
private:
    double alpha;
    double dx;
    
public:
    HeatEquationSystem(double _alpha, double _dx) : alpha(_alpha), dx(_dx) {}
    
    void operator()(const Eigen::VectorXd &u, Eigen::VectorXd &dudt, const double /* t */) 
    {
        const size_t size = u.size();

        //check size before any operations
        if (size < 3) {
            throw std::runtime_error("Vector size too small for finite difference computation.");
        }

        //finite differencing 
        for (size_t i = 1; i < size - 1; ++i) {
            dudt(i) = alpha * (u(i - 1) - 2.0 * u(i) + u(i + 1)) / (dx * dx);
        }

        //boundary conditions on the time derivative 
        dudt(0) = 0;
        dudt(size - 1) = 0;
    }
};

int main() 
{
    std::ofstream file("solution.txt"); 

    const double alpha = 0.1; 
    const size_t M = 10;       // number of spatial intervals
    const size_t N = 100;      // number of time steps
    const double dx = 1.0/M;   // spatial step size 
    const double dt = 0.001;   // time step size
    constexpr double pi = 3.14159265358979323846;
    const double t_start = 0; 
    
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(M+1, 0, 1); //spatial grid
    
    //solution vector
    Eigen::VectorXd u(M+1); 

    //initial conditions
    for (size_t i = 0; i <= M; ++i) {
        u(i) = sin(pi * x(i)); 
    }

    //boundary conditions
    u(0) = 0; 
    u(M) = 0; 

    //ODE system
    HeatEquationSystem system(alpha, dx); //odeint require a callable object 


    //setupp ODE solver
    odeint::runge_kutta4<Eigen::VectorXd, double, Eigen::VectorXd, double, odeint::vector_space_algebra> stepper;

    // Time integration
    double t = t_start;
    for (size_t n = 0; n < N; ++n) 
    {
        // Save current solution to file
        file << t << " ";
        for (size_t i = 0; i <= M; ++i) {
            file << u(i) << " ";
        }
        file << "\n";

        // Perform time step
        stepper.do_step(system, u, t, dt);
        t += dt;
    }

    file.close();
    return 0;
}