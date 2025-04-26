#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>

using namespace Eigen; 


double binomial_american_option(const double& T, const int& N, const double& sigma, const double& S0, const double& r, const bool& is_call, const double& K)
{
    auto start = std::chrono::high_resolution_clock::now(); 
    //discretize time 
    double dt = T /N; 

    //ups and downs factors 
    double u = exp(sigma * sqrt(dt)); 
    double d = 1.0/u; 

    //risk neautral probability 
    double p = (exp(r*dt)-d)/(u-d); 
 
    MatrixXd S = MatrixXd::Zero(N+1, N+1) ; //price tree 
    MatrixXd V = MatrixXd::Zero(N+1, N+1); //option value tree 

    //construct the binomial tree
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
    auto end = std::chrono::high_resolution_clock::now(); 
    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)   << "\n"; 


    return V(0,0); //option price at t=0

} 


double vector_binomial_american_option(const double& T, const int& N, const double& sigma, const double& S0, const double& r, const bool& is_call, const double& K)
{
    auto start = std::chrono::high_resolution_clock::now();

    using RowMatrix = Matrix<double, Dynamic, Dynamic, RowMajor>; //define a matrix with dimension defined at compile time that uses rowmajor storage (to improve cache locality)

    //discretize time 
    double dt = T /N; 

    //ups and downs factors 
    double u = exp(sigma * sqrt(dt)); 
    double d = 1.0/u; 

    //risk neautral probability 
    double p = (exp(r*dt)-d)/(u-d); 

    //discount factor (compute it once)
    double discount = exp(-r*dt); 


    RowMatrix S(N+1, N+1); //stock price tree 
    RowMatrix V(N+1, N+1); //option value tree 

    // auto idx = [] (int t, int j) {return t *(t+1)/2 + j; }; //assing lambda function to a varibale to make it a functor 

    //construct the stock price tree
    S(0,0) = S0; 
    for (int t=1; t<N+1; ++t)
    {
        S(t,0) = S(t-1, 0) *u; 
        S.block(t,1,1,t) = S.block(t-1, 0, 1, t).array() *d;  //vector operation 
    }

    //compute terminal payoff 
    if (is_call)
        V.row(N) = (S.row(N).array() - K).cwiseMax(0.0); //cwise is the element wise max 
    
    else    
        V.row(N) = (K - S.row(N).array()).cwiseMax(0.0); 


    //compute option value tree 
    for (int t=N-1; t>=0; t-- )
    {
        RowMatrix CV = (V.block(t+1, 0, 1, t+1)* p + V.block(t+1, 1, 1, t+1) *(1-p))*discount; 

        RowMatrix EV; 
        if(is_call)
            EV = (S.block(t,0,1,t+1).array()-K).cwiseMax(0.0); 

        else    
            EV = (K - S.block(t,0,1,t+1).array()).cwiseMax(0.0); 

        //compute option value at time t 
        V.block(t,0,1,t+1) = CV.array().max(EV.array());
    }
    auto end = std::chrono::high_resolution_clock::now(); 
    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)   << "\n"; 


    return V(0,0); 

}


int main ()
{
    //parameters 
    double S0 = 100.0;   //current stock price
    double K = 105.0;    //strike price
    double T = 1.0;      //time to expiration (1 year)
    double r = 0.05;     //risk-free interest rate (5%)
    double sigma = 0.2;  //volatility (20%)
    int N = 8000;         //number of time steps
    bool is_call = true; //true for call option, false for put option

    std::cout << "non vectorized binomial\n" <<  binomial_american_option(T, N, sigma, S0, r,  is_call, K) << "\n";
    std::cout << "vectorized binomial\n" << vector_binomial_american_option(T, N, sigma, S0, r,  is_call, K) << "\n"; 
}