#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <random>

using namespace Eigen; 

VectorXd gaussianNoise(const int& N)
{
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::normal_distribution<double> dis(0,1); 

    VectorXd noise (N);  

    for (int i =0; i< N; ++i)
        noise(i) = dis(gen); 

    return noise; 
}

double call_option_monteCarlo(const int& M, const double& S0, const double & K, const double& r , const double& sigma, const double& T)
{ 
    //compute stock price 
    double S_adjust = S0 * exp((r-0.5*sigma *sigma)*T) ; //stock price with only the drift 
    double payoff_sum =0.0; 
    double S_curr = 0.0; 
    auto noise = gaussianNoise(M); 

    for (int i =0; i< M ; ++i)
    {
        S_curr = S_adjust *exp(sigma*sqrt(T)*noise(i)); //stock price with the random shock 
        payoff_sum += std::max(S_curr-K, 0.0); 
    }

    return exp(-r*T)*payoff_sum/M; 

}

double put_option_monteCarlo(const int& M, const double& S0, const double & K, const double& r , const double& sigma, const double& T)
{ 
    //compute stock price 
    double S_adjust = S0 * exp((r-0.5*sigma *sigma)*T) ; //stock price with only the drift 
    double payoff_sum =0.0; 
    double S_curr = 0.0; 
    auto noise = gaussianNoise(M); 

    for (int i =0; i< M ; ++i)
    {
        S_curr = S_adjust *exp(sigma*sqrt(T)*noise(i)); //stock price with the random shock 
        payoff_sum += std::max(K-S_curr, 0.0); 
    }

    return exp(-r*T)*payoff_sum/M; 

}


int main() 
{
    int M = 10000000; //n simulations 
    double S = 100.0 ; //option price 
    double K = 100.0 ; //strike price 
    double r = 0.05 ; //risk free rate 
    double sigma = 0.2 ; //volatility of the stock  
    double T = 1.0 ; // expirity date 

    std::cout << "call option price: " << call_option_monteCarlo(M, S, K, r, sigma, T) << "\n"; 
    std::cout << "put option price: " << put_option_monteCarlo(M, S, K, r, sigma, T) << "\n"; 

}