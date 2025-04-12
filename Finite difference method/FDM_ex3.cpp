//solve PDE (Heat Map): ut = uxx 
//conditions:  initial value: u(x, 0) = sin(pi*x) , boundary value: go(t) = g1(t) =0 

#include <Eigen/Dense>
#include <iostream>
#include <fstream>

int main ()
{
    constexpr double pi = 3.14159265358979323846; 

    size_t M = 4; //space intervals 
    size_t N = 40;  // time intervals 

    double dx = 1.0/M; //space steps size 
    double dt = 1.0/N; //time steps size  
    std::cout << "dx: " << dx << " dt: " << dt << "\n"; 

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(M+1, 0,1); //space gridpoints 
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N+1, 0, 1); //time gridpoints

    Eigen::MatrixXd U (M+1, N+1);  //solution matrix (space, time)

    //intial conditions 
    for (int i =0; i< M+1; ++i)
        U(i, 0) = sin(pi * x(i)); 
        
    //boundary conditions 
    U.row(0) = Eigen::VectorXd::Zero(N+1); 
    U.row(M) = Eigen::VectorXd::Zero(N+1); 

    //check stability
    double stabilityConditions = dt / (dx *dx); 
    std::cout << "r: " << stabilityConditions << "\n"; 
    if (stabilityConditions > 0.5)
        std::cout << "Stability conditions not satisfied\n"; 

    
    //explicit Euler method
    double r = dt/ (dx*dx); 
    for (int j=0; j<N; ++j)
        for(int i =1; i< M; ++i)
            U(i, j+1) = U(i, j) +r*(U(i+1, j)-2*U(i, j) +U(i-1, j)); 

    //visualize
    std::cout << "solution\n" << U; 

    std::ofstream file("solution"); 

    if(file)
    {
        file << U; 
        file.close(); 
        std::cout << "\nsolution saved in the file\n"; 
    }

}
