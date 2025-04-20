#include <Eigen/Dense>
#include <iostream>
#include <fstream>

// PDE : ut = uxx
//conditions: initial values U(x, 0):  2x if 0 <= x <= 0.5 , 2(1-x) if 0.5 < x <= 2-2x
//boundary values: go(t) = g1(t) = 0

int main ()
{
    size_t M = 4; //space steps 
    size_t N = 20; //time steps 

    double dx = 1.0/M; //space steps size 
    double dt = 1.0/N; //time steps size 
    std::cout << "dx: " << dx << " dt: " << dt << "\n"; 

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(M+1, 0, 1); //space grid points 
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N+1, 0, 1); //time grid points

    Eigen::MatrixXd U (M+1, N+1); //solution matrix 

    //intial conditions 
    for (int i =0; i< M+1; ++i)
    {
        if (x(i)>= 0 && x(i)<= 0.5)
            U(i, 0) = 2*x(i); 
        
        else if (x(i) > 0.5 && x(i) <= (2-2*x(i)))
            U(i, 0) = 2*(1-x(i)); 
    }

    //boundary conditions 
    U.row(0) = Eigen::VectorXd::Zero(N+1);  
    U.row(M) = Eigen::VectorXd::Zero(N+1); 


    //check stability 
    double r = dt/(dx*dx); 
    if (r >= 0.5)
        std::cout << "Stability conditions not satisfied "; 

    //solve ODE with explicit Euler method
    for (int j=0; j< N; ++j)
        for (int i =1 ; i< M; ++i)
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