#include <Eigen/Dense>
#include <iostream>

//differential equation: u'' + 2u' - 3u = 9x  ,  [0,1]

//conditions: u'(0) = -4, u(1) = -2e^-3 + e - 5


int main()
{
    double a =0; 
    double b =1; 
    size_t N = 4; //n of steps
    double h = (b-a) /N; 
    double ub = exp(-3) +2*exp(1) -5; 
    std::cout << "h: " << h << "\n\n"; 

    double di = -2 - 3*h*h; 
    double wi = 1+h; 
    double vi = 1-h; 

    //discretize space 
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N+1, a, b); 

    //generate tridiagonal matrix A
    Eigen::VectorXd centraDiagonal = Eigen::VectorXd::Constant(N+1, di);  
    Eigen::VectorXd lowerDiagonal = Eigen::VectorXd::Constant(N, vi);  
    Eigen::VectorXd upperDiagonal = Eigen::VectorXd::Constant(N, wi);  

    Eigen::MatrixXd A(N+1, N+1); 
    A = centraDiagonal.asDiagonal();
    A.diagonal(-1) = lowerDiagonal; 
    A.diagonal(1) = upperDiagonal;  

    //generate B
    Eigen::VectorXd B = 9*h*h*x;

    //apply constrictions 
    A(0,0) = -2-3*h*h;
    A(0,1) = 2; 
    B(0) = (2*h -2*h*h)*-4;  
    
    A(N, N) = 1; 
    A(N, N-1) = 0; 
    B(N) = ub; 

    //visualize 
    std::cout << "A\n" << A << "\n\n"; 
    std::cout << "B\n" << B << "\n\n";

    //compute solution Ax = B
    Eigen::VectorXd solution = A.lu().solve(B); 
    std::cout << "solution\n" << solution << "\n\n"; 

}