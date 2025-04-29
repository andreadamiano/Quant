#include <iostream>
#include <cmath>

double norm_cdf(const int x)
{
    return 0.5 * (1.0 + erf(x/sqrt(2.0))); 
}

double di( const int i, const double S , const double K, const double r ,const double sigma , const double T)
{
    return (log(S/K) + (r+sigma*sigma/2.0)*T - pow(-1,i-1) )/ (sigma * sqrt(T)); 
}

double call_price(const double S , const double K, const double r , const double sigma , const double T)
{
    return S*norm_cdf(di(1, S, K, r, sigma, T)) - K*exp(-r*T) * norm_cdf(di(2, S, K, r, sigma, T)); 
}

double call_delta (const double S , const double K, const double r , const double sigma , const double T, const double delta_S )
{
    return (call_price(S + delta_S, K, r, sigma, T) - call_price(S -delta_S, K, r, sigma, T))/2.0*delta_S; 
}

double call_gamma (const double S , const double K, const double r , const double sigma , const double T, const double delta_S )
{
    return (call_price(S + delta_S, K, r, sigma, T) - 2*call_price(S , K, r, sigma, T) + call_price(S - delta_S, K, r, sigma, T)) /(delta_S*delta_S); 
}

int main ()
{
    double S = 100.0 ; // Option price
    double delta_S = 0.001 ; // Option price increment
    double K = 100.0 ; // Strike price
    double r = 0.05 ; // Risk−free rate (5%)
    double sigma = 0.2 ; // Volatility of the underlying (20%)
    double T = 1.0 ; // One year until expiry

    //compute greeeks 
    std::cout << "delta: " << call_delta(S, K, r, sigma, T, delta_S) << "\n"; 
    std::cout << "gamma: " << call_gamma(S, K, r, sigma, T, delta_S) << "\n"; 
}