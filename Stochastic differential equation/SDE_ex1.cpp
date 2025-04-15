#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

//SDE: dX(t) = −µX(t) dt + σ dWt

void EulerMayurama(const int& N, const double& T, Eigen::VectorXd X, const double& mean, const double& sigma, const double& X0)
{
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::normal_distribution<double> dis (0, 1);    

    double dW=0.0; 
    double dt = T/N; 

    X(0) = X0;
    for (int i =0; i< N; ++i)
    {
        //generate Brownian steps 
        dW= sqrt(dt) * dis(gen); 
        X(i+1) = X(i)-(mean*X(i))*dt + sigma *dW; 

    }
}

int main()
{
    const int N = 1000;
    const double T = 1.0;
    const double mean = 2.0;
    const double sigma = 0.3;
    const double X0 = 1.0;
    Eigen::VectorXd X (N+1); 
    std::ofstream file ("result.csv"); 

    EulerMayurama(N, T, X, mean, sigma, X0); 

    if(file)
    {
        file << X; 
    }

}