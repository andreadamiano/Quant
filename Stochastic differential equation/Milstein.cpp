#include <Eigen//Dense>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

//SDE: // SDE: dX_t = μ * X_t * dt + σ * X_t * dW_t (GBM)

using namespace Eigen; 

void Milstein(const double& T, const int& N, const double& mean, const double& sigma, MatrixXd& X , const double& X0)
{
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::normal_distribution<double> dis(0.0, 1.0); 
    double dt = T/N; 

    //Brownian Motion
    VectorXd dW(N); //Brownian increments
    VectorXd W(N+1); //cumulative BM 
    W(0) =0; 

    //initialize solution
    X(0,0) = X0; //Millstein solution
    X(0,1) = X0; //analytical solution 

    for (int i =0; i< N; i++)
    {
        dW(i) = dis(gen) * sqrt(dt); //brownian increment 
        W(i+1) = W(i) + dW(i); 

        //numerical solution
        X(i+1, 0) = X(i,0) + mean*X(i,0)*dt + sigma*X(i, 0)* dW(i) 
        + 0.5*sigma*sigma*X(i) * (dW(i)*dW(i) -dt); //Millstein correction term 

        //analytical solution
        X(i+1, 1) = X0* exp((mean-0.5*sigma*sigma)*(i+1)*dt + sigma* W(i+1)); 
    }

}


int main ()
{
    // Parameters
    const double T = 1.0; //time horizon
    const int N = 1000;  //time steps
    const double X0 = 1.0; //initial condition
    const double mean = 0.1; //drift
    const double sigma = 0.2;  //volatility 

    MatrixXd X (N+1,2);
    std::ofstream file ("result.csv"); 
    
    Milstein(T, N, mean, sigma, X, X0); 

    if (file)
    {
        file << X; 
    }


}