#include <Eigen/Dense>
#include <fstream>
#include <iostream>

//Black and Scholes PDE: ∂V/∂t + (1/2)σ²S²(∂²V/∂S²) + rS(∂V/∂S) - rV = 0

// Where:
// - V(S,t) is the option price as function of underlying price S and time t
// - σ is the volatility of the underlying asset
// - r is the risk-free interest rate
// - ∂V/∂t is the time derivative of the option price
// - ∂V/∂S is the first derivative with respect to the underlying price
// - ∂²V/∂S² is the second derivative with respect to the underlying price

//conditions 
//initial condition: Vi(0) = max(si - K, 0)
//boundary conditions: v0(t) = 0 , Vk(t) = 200−100e^(-0.05*(1-t))


int main ()
{
    double T = 1.0; //1 year
    double S_max = 200.0; //max stock price  
    double sigma = 0.3; //volatility
    double r = 0.05; //risk free interest rate 
    double K =100; //strike price 
    size_t M = 100; //space steps 
    size_t N = 1000; //time steps

    double dS = S_max/M;
    double dt = T/N; 
     
    std::cout << "dx: " << dS << " dt: " << dt; 

    Eigen::VectorXd S = Eigen::VectorXd::LinSpaced(M+1, 0, S_max); //space grid points 
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N+1, 0, T);  //time grid points 

    Eigen::MatrixXd U (M+1, N+1); //solution matrix (space, time)

    //conditions
    //intial condition
    for (int i =0; i< M+1; ++i)
    {
        U(i, 0) = std::max(S(i)-K, 0.0); 
    }

    //boundary condition
    for (int j=0; j< N+1; ++j)
        U(M, j) = S_max - K * exp(-r*(T-t(j))); //boundary condition
    U.row(0) = Eigen::VectorXd::Zero(N+1); 


    //define coefficients of the ODE 
    Eigen::VectorXd alpha (M+1); 
    Eigen::VectorXd beta (M+1); 
    Eigen::VectorXd gamma (M+1); 

    for (int i=1; i<M; ++i) //computed for interior points 
    {
        alpha(i) = -0.5 * sigma *sigma * S(i) * S(i) * 1/(dS*dS) + (r*S(i))/(2*dS); 
        beta(i) = sigma *sigma * S(i) * S(i) * 1/(dS*dS) + r;   
        gamma(i) = -0.5 *sigma *sigma * S(i) * S(i) * 1/(dS*dS) - (r*S(i))/(2*dS); 
    }


    //check stability
    if (dt > (dS*dS)/(sigma*sigma*S_max*S_max+r*dS*dS))
        std::cout << "Stability conditions not satisfied "; 

    
    //solve ODE with explicit Euler method
    for (int j=0; j< N; ++j)
        for (int i =1 ; i< M; ++i)
        U(i, j+1) = U(i, j) + dt * (alpha(i)*U(i-1, j) + beta(i)* U(i, j)+ gamma(i)*U(i+1, j)); 


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