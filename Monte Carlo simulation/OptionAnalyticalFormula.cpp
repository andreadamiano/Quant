#include <iostream>
#include <cmath>

//Black and Scholes analytical formula

constexpr double gamma = 0.2316419; 
constexpr double b1 = 0.319381530 ; 
constexpr double b2 = -0.356563782; 
constexpr double b3 = 1.781477937; 
constexpr double b4 = -1.821255978; 
constexpr double b5 = 1.330274429 ; 

double norm_pdf(const double& x)
{
    return (1.0/(pow(2*M_PI, 0.5))* exp(-0.5*(x*x))); 
}

double norm_cdf (const double& x)
{
    double k =  1/(1+gamma*x);
    double k_sum = k * (b1 + k * (b2 + k * (b3 + k * (b4 + k * b5))));

    if(x>=0)
        return 1-norm_cdf(x)*k_sum; 

    else    
        return 1-norm_cdf(-x); 


}

int main ()
{

}