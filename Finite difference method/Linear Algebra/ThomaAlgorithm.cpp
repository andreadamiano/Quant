#include <Eigen/Dense>
#include <iostream>

using namespace Eigen; 

void ThomasAlgo(const VectorXd& a,const VectorXd& b, const VectorXd& c, const VectorXd& d, VectorXd& f)
{
    size_t size = d.size(); //size of the matrix 

    VectorXd c_star = VectorXd::Zero(size); 
    VectorXd d_star = VectorXd::Zero(size); 

    //initialize the beginning elements of c_star and d_star
    c_star(0) = c(0)/b(0); 
    d_star(0) = d(0)/b(0); 

    //forward step 
    for (int i =1; i < size; ++i )
    {
        double m = 1/ (b(i)-c_star(i-1)*a(i)); //elimination factor 
        
        c_star(i) = c(i) *m;  
        d_star(i) = (d(i) -d_star(i-1)*a(i)) * m; 
    }

    //backward substitution
    f(size-1) = d_star(size-1);
    for (int i = size-2; i >= 0; --i)
        f(i) = d_star(i) - c_star(i)*f(i+1);

}


int main ()
{
    //use thomas algorithm to solve tridiagonal system

    size_t N = 13; //space steps
    double dx = 1.0/N; //space steps size 
    double dt = 0.001; //time step size 
    double r = dt /(dx*dx); 


    //vector of coefficients 
    Eigen::VectorXd a = Eigen::VectorXd::Constant(N, -r/2.0); 
    Eigen::VectorXd b = Eigen::VectorXd::Constant(N, 1.0 + r); 
    Eigen::VectorXd c = Eigen::VectorXd::Constant(N, -r/2.0); 
    Eigen::VectorXd d = Eigen::VectorXd::Constant(N, 0.0); 
    Eigen::VectorXd f = Eigen::VectorXd::Constant(N, 0.0); 

    //intial value
    f(5) = 1; 
    f(6)= 2; 
    f(7) = 1; 

    //fill vector d 
    for (int i =1; i<N-1; ++i)
        d(i) = r * 0.5 *f(i+1) + (1.0-r)*f(i) +r *0.5*f(i-1); 
    
    //solve tridiagonal system 
    ThomasAlgo(a,b,c,d,f); 

    std::cout << "solution: \n" << f.transpose(); 

}

