#include "CallOption.h"
#include "Black&scholes.h"
#include <memory>
#include <stdexcept>
#include <cmath>

double norm_cdf(const double x)
{
    return 0.5 * (1.0+ erf(x/sqrt(2.0))); 
}

//implementation details 

CallOption::CallOption(double strike, double maturity): strike(strike), maturity(maturity) {}

double CallOption::price(const BlackScholesModel& model) 
{
    //price the option using black and scholes analytical solution 
    double S = model.stockPrice;    
    double K = this->strike; 
    double T = this->maturity - model.date; //time to expiration
    double r = model.riskfreeRate; 
    double sigma = model.volatility; 

    if (T < 0.0) //if time is expired 
        return 0.0; 

    //black and scholes formula 
    double d1 = (log(S/K) + (r+ sigma*sigma*0.5)*T)/(sigma*sqrt(T)); 
    double d2 = d1 - sigma*sqrt(T); 

    return S * norm_cdf(d1) - K *exp(-r*T) * norm_cdf(d2); 
}

std::unique_ptr<CallOption> CallOption::newInstance(const double strike, const double maturity)
{
    if (strike <=0 || maturity <=0)
        throw std::runtime_error("error, invalid strike/maturity"); 

    return std::unique_ptr<CallOption> (new CallOption(strike, maturity)); //return a unique pointer to a new object 
}