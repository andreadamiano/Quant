#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen; 


double binomial_american_option(const double& T, const int& N, const double& sigma, const double& S0, const double& r, const bool& is_call, const double& K)
{
    //discretize time 
    double dt = T /N; 

    //ups and downs factors 
    double u = exp(sigma * sqrt(dt)); 
    double d = 1.0/u; 

    //risk neautral probability 
    double p = (exp(r*dt)-d)/(u-d); 


    //construct the binomial tree 
    MatrixXd S = MatrixXd::Zero(N+1, N+1) ; //price tree 
    MatrixXd V = MatrixXd::Zero(N+1, N+1); //option value tree 

    S(0,0) = S0; //initial stock price 
    for (int t=1; t <N+1; ++t)
    {
        S(t, 0) = S(t-1, 0) * u; 
        for (int j=1; j<=t; ++j)
        {
            S(t, j) = S(t-1, j-1) * d; 
        }
    }//result to a lower triangular matrix 

    //compute payoff at maturity for each leaf 
    for (int j=0; j<N+1; j++)
    {
        if(is_call)
            V(N, j) = std::max(S(N, j)-K, 0.0); 

        else 
            V(N, j) = std::max(K-S(N, j), 0.0); 
    }

    //backward induction of the option price , with early exercize 
    for (int t = N-1; t>=0; --t)
    {
        for (int j=0; j<=t; ++j)
        {
            double CV = exp(-r*dt)*(V(t+1,j)*p + (1-p)*V(t+1,j+1)); //continuation value (if option is not exercised)

            double EV = 0.0; 
            if  (is_call)
                EV = std::max(S(t,j)- K, 0.0); 
            else    
                EV = std::max(K - S(t,j), 0.0); 

            V(t,j) = std::max(CV, EV); 
        }
    }

    return V(0,0); //option price at t=0

} 


int main ()
{
    //parameters 
    double S0 = 100.0;   //current stock price
    double K = 105.0;    //strike price
    double T = 1.0;      //time to expiration (1 year)
    double r = 0.05;     //risk-free interest rate (5%)
    double sigma = 0.2;  //volatility (20%)
    int N = 4000;         //number of time steps
    bool is_call = true; //true for call option, false for put option

    std::cout << binomial_american_option(T, N, sigma, S0, r,  is_call, K);
}