#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace Eigen;

void ThomasAlgo(const VectorXd& a, const VectorXd& b, const VectorXd& c, const VectorXd& d, VectorXd& f) {
    size_t size = d.size();
    VectorXd c_star(size);
    VectorXd d_star(size);

    // Forward sweep
    c_star(0) = c(0)/b(0);
    d_star(0) = d(0)/b(0);
    
    for (int i = 1; i < size; ++i) 
    {
        double m = 1.0/(b(i) - a(i)*c_star(i-1));
        c_star(i) = c(i)*m;
        d_star(i) = (d(i) - a(i)*d_star(i-1))*m;
    }

    // Backward substitution
    f(size-1) = d_star(size-1);
    for (int i = size-2; i >= 0; --i) {
        f(i) = d_star(i) - c_star(i)*f(i+1);
    }
}

int main() {
    std::ofstream file("solution.csv");
    
    size_t M = 50;         // Space steps
    double dx = 1.0/M;      // Space step size
    double dt = 0.01;      // Time step size
    double c = 2.0;         // Wave speed
    size_t N = 2000;      // Time steps
    const double pi = 3.14159265358979323846;

    VectorXd x = VectorXd::LinSpaced(M+1, 0, 1);

    // Initialize solution matrix
    MatrixXd U = MatrixXd::Zero(M+1, N+1);

    // Initial conditions
    for (int i = 0; i < M+1; ++i) 
    {
        U(i, 0) = sin(pi*x(i));
    }

    // Initial velocity condition (first time step)
    double r = c*dt/dx;
    for (int i = 1; i < M; ++i) 
    {
        U(i, 1) = U(i, 0) + (r*r/2.0)*(U(i+1, 0) - 2*U(i, 0) + U(i-1, 0));
    }

    // Boundary conditions
    for (int n = 0; n < N+1; ++n) 
    {
        U(0, n) = 0.0;
        U(M, n) = 0.0;
    }

    // Crank-Nicolson coefficients
    double alpha = r*r/4.0;
    VectorXd a = VectorXd::Constant(M-1, -alpha);    //lower diagonal
    VectorXd b = VectorXd::Constant(M-1, 1.0+2.0*alpha); //main diagonal
    VectorXd c_vec = VectorXd::Constant(M-1, -alpha);    //upper diagonal
    VectorXd d_vec = VectorXd::Zero(M-1);           //right-hand side
    VectorXd u_sol = VectorXd::Zero(M-1);           //solution vector

    //time-stepping loop (for internal points only)
    for (int n = 1; n < N; ++n) 
    {
        for (int i = 1; i < M; ++i) 
        {
            d_vec(i-1) = alpha*U(i-1, n) + (1.0-2.0*alpha)*U(i, n) + alpha*U(i+1, n); //construct right-hand side
        }

        //solve tridiagonal system
        ThomasAlgo(a, b, c_vec, d_vec, u_sol);

        // Update solution
        for (int i = 1; i < M; ++i) 
        {
            U(i, n+1) = u_sol(i-1);
        }
    }

    //visualize
    for (int n = 0; n <= N; n++) 
    {
        for (int i = 0; i <= M; i++) 
        {
            file << U(i, n) << " ";
        }
        file << "\n";
    }
    file.close();

    std::cout << "Simulation completed. Results saved to solution.csv" << std::endl;
    return 0;
}