#include "CallOption.h"
#include <iostream>

int main()
{
    //define black and scholes model 
    auto model = BlackScholesModel(); 
    model.stockPrice = 150.0; 
    model.volatility = 0.25; 
    model.riskfreeRate = 0.83; 
    model.date = 0.0; 

    double strike = 100.0; 
    double maturity = 1.0; 
    auto option = CallOption::newInstance(strike, maturity); 
    std::cout << option->price(model); 
}