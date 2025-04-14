//integral between 0 and 1 of x^2

#include <iostream>
#include <cmath>
#include <random>

double fun(const double& x)
{
    return x*x; 
}


int main ()
{
    constexpr double N =1000.0; 
    constexpr double a =0.0; 
    constexpr double b =1.0; 


    std::random_device rd ; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution dis(a,b); 

    double sum(0); 
    for (int i=0; i < N; i++)
    {
        sum += fun(dis(gen));    
    }

    double integral = (b-a)*(sum/N); //solution 0.333

    std::cout << integral; 
}