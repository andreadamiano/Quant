#include <Eigen/Dense>
#include <iostream>
#include <cmath>

//disccretise ODE
//u'' + 2u' - 3u = 9x on the interval [0,1]

//boundaary conditions: u(0) = 1 ,  u(1) = e^-3 +2e - 5 circa 0.486351




int main()
{
    double a =0; 
    double b =1; 

    double ua = 1.0; 
    double ub = exp(-3) + 2*exp(1)-5; 

    size_t N = 4; 
    double h = (b-a)/N; 
    
    std::cout << "h: " << h << "\n"; 

    double di = -2-3*h*h; 
    double vi = 1-h; 
    double wi = 1+h; 

    //discretize the space 
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N+1, a, b); 

    //generate tridiagonal matrix A 
    Eigen::VectorXd centralDiagonal =  Eigen::VectorXd::Constant(N+1, di);
    Eigen::VectorXd lowerDiagonal =  Eigen::VectorXd::Constant(N, vi);
    Eigen::VectorXd upperDiagonal =  Eigen::VectorXd::Constant(N, wi);

    Eigen::MatrixXd A(N+1, N+1); 
    A = centralDiagonal.asDiagonal(); 
    A.diagonal(1) = upperDiagonal; 
    A.diagonal(-1) = lowerDiagonal; 

    //generate b
    Eigen::VectorXd B = 9*x*h*h; 

    //modify first equation to impose the left boundary cosntriction
    A(0,0) = 1; 
    A(0,1) =0; 
    B(0) = ua; 

    //modify second equation to impose right boundary cosntriction
    A(N, N) =1; 
    A(N, N-1)=0; 
    B(N) = ub; 

    std::cout << "A\n" << A << "\n\n"; 
    std::cout << "b\n" << B << "\n\n"; 


    //solve Ax = B
    Eigen::VectorXd solution = A.lu().solve(B); 
    std::cout << "solution\n" << solution << "\n\n";

}