#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <random>
#include <omp.h>
#include <chrono>

using namespace Eigen; 

VectorXd gaussianNoise(const int& N)
{

    VectorXd noise (N);  

    #pragma omp parallel num_threads(4)
    {
        std::random_device rd; 
        std::mt19937 gen(rd()); 
        std::normal_distribution<double> dis(0,1); 

        #pragma omp for
        for (int i =0; i< N; ++i)
            noise(i) = dis(gen); 
    }

    return noise; 
}

void parallel_monte_carlo_call_option (const int num_sims, const double S , const double K, const double r,const double sigma , const double T, 
const double delta_S, double& priceSp , double& priceS , double& priceSm)
{
    //generate 3 different stock paths for each increment / decrement 
    const double Sp_adjust = (S + delta_S) * exp(T*(r-0.5*sigma*sigma)); 
    const double S_adjust = (S) * exp(T*(r-0.5*sigma*sigma)); 
    const double Sm_adjust = (S - delta_S) * exp(T*(r-0.5*sigma*sigma)); 

    double payoff_sum_p = 0.0;
    double payoff_sum = 0.0;
    double payoff_sum_m = 0.0;

    const double vol_factor = sigma * sqrt(T);
    const int batch_size = 32768;  

    const VectorXd noise = gaussianNoise(num_sims);

    omp_set_nested(1); // enable nested parallelism 
    #pragma omp parallel sections num_threads(3)
    {
        //Sp price path 
        #pragma omp section
        {
            #pragma omp parallel num_threads(4) reduction(+:payoff_sum_p)
            {
                VectorXd S_curr(batch_size); //one vector per thread 
                VectorXd payoff(batch_size);

                #pragma omp for schedule(static)
                for (int batch_start = 0; batch_start < num_sims; batch_start += batch_size)
                {
                    const int current_batch_size = std::min(batch_size, num_sims - batch_start);
                    auto current_noise = noise.segment(batch_start, current_batch_size);
                    
                    S_curr.head(current_batch_size) = Sp_adjust * (vol_factor * current_noise).array().exp();
                    payoff.head(current_batch_size) = (S_curr.head(current_batch_size).array() - K).cwiseMax(0.0);
                    payoff_sum_p += payoff.head(current_batch_size).sum();
                }
            }
        }

        //S price path 
        #pragma omp section
        {
            #pragma omp parallel num_threads(4) reduction(+:payoff_sum)
            {
                VectorXd S_curr(batch_size); //one vector per thread 
                VectorXd payoff(batch_size);

                #pragma omp for schedule(static)
                for (int batch_start = 0; batch_start < num_sims; batch_start += batch_size)
                {
                    const int current_batch_size = std::min(batch_size, num_sims - batch_start);
                    auto current_noise = noise.segment(batch_start, current_batch_size);
                    
                    S_curr.head(current_batch_size) = S_adjust * (vol_factor * current_noise).array().exp();
                    payoff.head(current_batch_size) = (S_curr.head(current_batch_size).array() - K).cwiseMax(0.0);
                    payoff_sum += payoff.head(current_batch_size).sum();
                }
            }
        }

        //Sm price path 
        #pragma omp section
        {
            #pragma omp parallel num_threads(4) reduction(+:payoff_sum_m)
            {
                VectorXd S_curr(batch_size); //one vector per thread 
                VectorXd payoff(batch_size);

                #pragma omp for schedule(static)
                for (int batch_start = 0; batch_start < num_sims; batch_start += batch_size)
                {
                    const int current_batch_size = std::min(batch_size, num_sims - batch_start);
                    auto current_noise = noise.segment(batch_start, current_batch_size);
                    
                    S_curr.head(current_batch_size) = Sm_adjust * (vol_factor * current_noise).array().exp();
                    payoff.head(current_batch_size) = (S_curr.head(current_batch_size).array() - K).cwiseMax(0.0);
                    payoff_sum_m += payoff.head(current_batch_size).sum();
                }
            }
        }
    }

    priceSp = payoff_sum_p / static_cast<double> (num_sims) * exp(-r*T); 
    priceS = payoff_sum / static_cast<double> (num_sims) * exp(-r*T); 
    priceSm = payoff_sum_m / static_cast<double> (num_sims) * exp(-r*T); 

}

void vectorized_monte_carlo_call_price(const int num_sims, const double S, const double K, const double r, const double sigma, const double T,
const double delta_S, double& priceSp, double& priceS, double& priceSm)
{
    //pre compute common terms
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol_factor = sigma * sqrt(T);
    
    //initial prices
    const double Sp_adjust = (S + delta_S) * exp(drift); 
    const double S_adjust = S * exp(drift);
    const double Sm_adjust = (S - delta_S) * exp(drift);
    
    const VectorXd noise = gaussianNoise(num_sims);
    const VectorXd exp_noise = (vol_factor * noise).array().exp();

    //terminal prices for the three cases
    VectorXd Sp_curr = Sp_adjust * exp_noise;
    VectorXd S_curr = S_adjust * exp_noise;
    VectorXd Sm_curr = Sm_adjust * exp_noise;

    //payoffs (call option payoffs)
    VectorXd payoffs_p = (Sp_curr.array() - K).cwiseMax(0.0);
    VectorXd payoffs = (S_curr.array() - K).cwiseMax(0.0);
    VectorXd payoffs_m = (Sm_curr.array() - K).cwiseMax(0.0);

    //discounted expected payoffs
    priceSp = payoffs_p.mean() * exp(-r * T);
    priceS = payoffs.mean() * exp(-r * T);
    priceSm = payoffs_m.mean() * exp(-r * T);
}

int main ()
{
    double S = 100.0 ; // Option price
    double delta_S = 0.001 ; // Option price increment
    double K = 100.0 ; // Strike price
    double r = 0.05 ; // Risk−free rate (5%)
    double sigma = 0.2 ; // Volatility of the underlying (20%)
    double T = 1.0 ; // One year until expiry
    int num_sims = 100000000; 

    double priceS =0.0; 
    double priceSp =0.0; 
    double priceSm =0.0; 

    //nested parallel loops 
    auto start = std::chrono::high_resolution_clock::now(); 
    parallel_monte_carlo_call_option(num_sims, S, K, r, sigma, T, delta_S, priceSp, priceS, priceSm); 
    auto end = std::chrono::high_resolution_clock::now(); 

    std::cout << "priceSp: " << priceSp << " priceS: " << priceS << " priceSm: " << priceSm << "\n";
    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds> (end-start) << "\n\n";  


    //vectorized operations
    start = std::chrono::high_resolution_clock::now(); 
    vectorized_monte_carlo_call_price(num_sims, S, K, r, sigma, T, delta_S, priceSp, priceS, priceSm); 
    end = std::chrono::high_resolution_clock::now(); 

    std::cout << "priceSp: " << priceSp << " priceS: " << priceS << " priceSm: " << priceSm << "\n";
    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds> (end-start) << "\n";  

}