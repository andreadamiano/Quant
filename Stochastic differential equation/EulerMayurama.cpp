#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

//simulate Heston Model 

void EulerMayurana(const double& rho, const double& T, const int& N,  Eigen::VectorXd& S,  Eigen::VectorXd& v , const double& mean, const double & theta , const double& k, const double & xi, const double& S0, const double& v0) 
{
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::normal_distribution<double> dis (0,1);
    double dt = T/N; //time step size

    S(0) = S0; 
    v(0) =v0; 


    for (int i =0; i< N; ++i)
    {
        //introduce correlation with Cholesky decomposition 
        double Z1 = dis(gen); 
        double Z2 = rho *Z1 + sqrt(1-rho*rho) * dis(gen); 

        //correlated Brownian increments
        double dW1 = sqrt(dt) * Z1;
        double dW2 = sqrt(dt)*Z2; 

        //update volatility
        v(i+1) = v(i) + k*(theta - v(i)) * dt + xi * sqrt(std::max(v(i), 0.0)) * dW2;

        //prevent negative volatility 
        v(i+1) = std::max(v(i+1), 0.0);
        
        //update stock price 
        S(i+1) = S(i) + mean*S(i)*dt+sqrt(v(i))*S(i)*dW1; 
    }



}

int main()
{
    double S0 = 100.0; //initial stock price 
    double v0 = 0.04; //volatility
    double mean = 0.0; //drift
    double k = 1.0; //mean reversion speed
    double theta = 0.04; //long-term
    double xi = 0.2; //vol-of-vol
    double rho = -0.7; //correlation
    double T = 1.0;//time horizon
    int N = 252; //n of steps 

    std::ofstream file ("result.csv"); 

    Eigen::VectorXd S(N+1); 
    Eigen::VectorXd v(N+1); 

    EulerMayurana(rho, T, N, S, v, mean, theta, k, xi, S0, v0); 

    if (file)
    {
        for (int i=0; i < N+1; ++i)
        {
            file << S(i) << " " << v(i) << "\n"; 
        }
    }

    file.close(); 
}