#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <fstream>

using namespace Eigen; 

VectorXd NormalNoise(const int& N)
{
    std::random_device rd; 
    std::mt19937 gen(rd()); //generate uniform distributed random noise 
    std::normal_distribution<double> dis(0, 1); //mean , std 

    VectorXd noise (N); 
    for (int i =0; i<N; ++i)
        noise(i) = dis(gen); 

    return noise; 
}

//it takes the time of the simulation and N of time steps  
VectorXd GeometricBrownianMotion (const double& T, const int& N, const double& S0, const double& mean, const double& std )
{
    double dt = T/N; 
    VectorXd W(N+1); //geometric brownian motion
    VectorXd B(N+1); //generic brownian motion
    auto dW = NormalNoise(N); 

    //initialize brownian motion and GBM
    B(0) =0; 
    W(0) =S0; 

    for (int i =1; i< N+1; ++i)
    {
        B(i) = B(i-1)+(mean - 0.5 * std * std)*dt +std*sqrt(dt)*dW(i-1); 
        W(i) = S0 *exp(B(i)); 
    }

    return W; 

}

int main ()
{
    //parameters 
    double T = 1.0; 
    int N = 100; //time steps 
    double S0 = 100; //intial price
    double mean = 0.0 ; //drift 
    double std = 0.2;  //volatility
    int M = 30; //n simulations 


    MatrixXd result (M, N+1); 
    std::ofstream file ("GBM.csv"); 
    
    for (int i =0; i<M; ++i)
        result.row(i) = GeometricBrownianMotion(T, N, S0, mean, std); 

    if (file)
    {   
        file << result; 
    }
}   