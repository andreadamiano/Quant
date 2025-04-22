#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <fstream>

//generates 2 uniform random points 
// double boxMullerTransform()
// {
//     static std::random_device rd ; 
//     static std::mt19937 gen(rd()); 
//     static std::uniform_real_distribution dis(0.0, 1.0); 

//     double u1 = dis(gen); 
//     double u2 = dis (gen); 

//     double z1 = sqrt(-2*u1)*cos(u2); 
//     double z2 = sqrt(-2*u1)*sin(u2); 
// }

Eigen::VectorXd gaussianNoise(const int& N, const double& mean, const double& std)
{
    static std::random_device rd ; 
    static std::mt19937 gen(rd()); 
    static std::normal_distribution<double> dis(mean, std); //mean , std 

    Eigen::VectorXd noise (N); 
    for (int i =0; i< N; ++i)   
        noise(i) = dis(gen); 

    return noise; 

}

Eigen::VectorXd brownianMotion(const double& T , const int& N)
{
    double dt = T/N; 
    Eigen::VectorXd W = Eigen::VectorXd::Zero(N+1); 

    auto dW = gaussianNoise(N, 0, 1); 

    for (int i =1; i<N+1; ++i)
        W(i) = W(i-1) + sqrt(dt)*dW(i-1); 


    return W; 

}

int main ()
{
    double T = 1.0;      //total time
    int N = 1000;        //time steps 
    int M = 10;  //n of paths 
    std::ofstream file ("BM.csv"); 
    Eigen::MatrixXd matrix (M, N+1); 

    for (int i = 0; i< M; ++i)
        matrix.row(i) = brownianMotion(T, N);

    if (file )
        file << matrix; 

    return 0;
}   