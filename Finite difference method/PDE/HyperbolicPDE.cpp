#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>


//PDE: utt = c^2 * uxx

//boundary conditions: u(0, t) = 0 , u(1, t) =0
//initial conditions: u(x, 0) =sin (pi*x)  , ut(x, 0) =0

using namespace Eigen; 

#define pi 3.14159265358979323846

MatrixXd Upwind (const int& N, const int& M, const double& c)
{
    MatrixXd U = MatrixXd::Zero(M+1, N+1);  //time × Space
    MatrixXd V = MatrixXd::Zero(M+1, N+1);
    MatrixXd W = MatrixXd::Zero(M+1, N+1);

    double dt = 1.0 /M; //time steps 
    double dx = 1.0 / N; //space steps 
    double CFL = c * dt / dx;  
    
    if (CFL >= 1.0)
    {
        std::cout << "Instability"; 
        return MatrixXd(); 
    }

    //initial displacement 
    for (int j = 0; j < N+1; ++j) 
    {
        U(0, j) = 100 * sin(pi *j * dx);  
    }

    //initial velocity (already zero)

    //initial stress w(x,0) = c* ux(x, 0)
    for (int j = 1; j < N; ++j) 
    {
        W(0, j) = c * (U(0, j+1) - U(0, j-1)) / (2 * dx); //central difference 
    } 

    //one sided difference for boundary points 
    W(0, 0) = c * (U(0, 1) - U(0, 0)) / dx;      
    W(0, N) = c * (U(0, N) - U(0, N-1)) / dx;


    //time-stepping
    for (int i = 0; i < M; ++i) 
    {
        //update velocity (except boundaries)
        for (int j = 1; j < N; ++j) {
            V(i+1, j) = V(i, j) + dt * (c * (W(i, j+1) - W(i, j-1)) / (2 * dx));
        }
        
        //update stress (except boundaries)
        for (int j = 1; j < N; ++j) {
            W(i+1, j) = W(i, j) + dt * (c * (V(i+1, j+1) - V(i+1, j-1)) / (2 * dx));
        }
        
        //boundary conditions for stress
        W(i+1, 0) = c * (U(i, 1) - U(i, 0)) / dx;
        W(i+1, N) = c * (U(i, N) - U(i, N-1)) / dx;
        
        //update displacement
        for (int j = 0; j <= N; ++j) {
            U(i+1, j) = U(i, j) + dt * V(i+1, j);
        }
        
        // Boundary conditions for displacement
        U(i+1, 0) = 0;
        U(i+1, N) = 0;
    }

    return U; 
}

//PDE : utt = c^2 * (uxx + uyy)
//initial conditions: u(x,y,0) = sin(πx)*sin(πy) , ut(x,y,0) =0
//boundary conditions: u(0,y,t) = u(1,y,t) = u(x,0,t) = u(x,1,t) = 0

std::vector<MatrixXd> Leapfrog(const int& Nx, const int& Ny, const int& Nt, const double& c)
{
    //discretize time and space 
    double dx = 1.0/Nx; 
    double dy = 1.0/Ny; 
    double dt = 1.0 /Nt; 
    double CFL = c * dt * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy));

    std::vector<MatrixXd> U(Nt+1, MatrixXd::Zero(Ny+1, Nx+1));  

    //check stability conditions 
    if (CFL > 1.0)
    {
        std::cout << "Instability"; 
        return {}; //return empty vector 
    }

    //initial conditions 
    for (int i =1; i< Ny; ++i)
    {
        for (int j =1; j <Nx; ++j)
        {
            U[0](i,j) = sin (pi*i*dy) * sin (pi*j*dx); 
        }
    }

    //first step 
    for (int i =1; i< Ny; ++i)
    {
        for (int j =1; j <Nx; ++j)
        {
            // double laplacian = (U[0](i+1, j)-2.0*U[0](i, j)+U[0](i-1,j))/(dy*dy) + (U[0](i, j+1)-2.0*U[0](i, j)+U[0](i,j-1))/(dx*dx) ; 
            // U[1](i, j) = U[0](i, j) + 0.5 * c*c*dt*dt * laplacian; 

            U[1](i,j) = sin(pi * j * dx) * sin(pi * i * dy) * cos(sqrt(2)* pi * c * dt);
        }
    }

    //boundary conditions (already taken into account by setting all to zero)

    //time loop (for internal points)
    // for (int n =1 ; n < Nt; n++)
    // {
    //     for (int i =1; i< Ny; ++i)
    //     {
    //         for (int j =1; j <Nx; ++j)
    //         {   
    //             double laplacian = (U[n](i+1,j)-2*U[n](i,j)+U[n](i-1,j))/(dx*dx) + (U[n](i,j+1)-2*U[n](i,j)+U[n](i,j-1))/(dy*dy);
    //             U[n+1](i,j) = 2*U[n](i,j) - U[n-1](i,j) + c*c*dt*dt*laplacian;
    //         }
    //     }

    //     //boundary conditions
    //     for (int i = 0; i <= Ny; ++i) 
    //     {
    //         U[n+1](i,0) = 0;      // x=0 boundary
    //         U[n+1](i,Nx) = 0;     // x=1 boundary
    //     }

    //     for (int j = 0; j <= Nx; ++j) 
    //     {
    //         U[n+1](0,j) = 0;      // y=0 boundary
    //         U[n+1](Ny,j) = 0;     // y=1 boundary
    //     }
    // }


    // Time loop
    for (int n = 1; n < Nt; ++n) {
        // Inner points (4th-order)
        for (int i = 2; i < Ny - 1; ++i) {
            for (int j = 2; j < Nx - 1; ++j) {
                double laplacian_x = (-U[n](i + 2, j) + 16 * U[n](i + 1, j) - 30 * U[n](i, j) + 16 * U[n](i - 1, j) - U[n](i - 2, j)) / (12 * dx * dx);
                double laplacian_y = (-U[n](i, j + 2) + 16 * U[n](i, j + 1) - 30 * U[n](i, j) + 16 * U[n](i, j - 1) - U[n](i, j - 2)) / (12 * dy * dy);
                U[n + 1](i, j) = 2 * U[n](i, j) - U[n - 1](i, j) + c * c * dt * dt * (laplacian_x + laplacian_y);
            }
        }

        // Near-boundary points (2nd-order)
        for (int i = 1; i < Ny; ++i) {
            for (int j = 1; j < Nx; ++j) {
                if (i < 2 || i > Ny - 2 || j < 2 || j > Nx - 2) {
                    double laplacian_x = (U[n](i + 1, j) - 2 * U[n](i, j) + U[n](i - 1, j)) / (dx * dx);
                    double laplacian_y = (U[n](i, j + 1) - 2 * U[n](i, j) + U[n](i, j - 1)) / (dy * dy);
                    U[n + 1](i, j) = 2 * U[n](i, j) - U[n - 1](i, j) + c * c * dt * dt * (laplacian_x + laplacian_y);
                }
            }
        }

        // Boundary conditions
        for (int i = 0; i <= Ny; ++i) {
            U[n + 1](i, 0) = 0;
            U[n + 1](i, Nx) = 0;
        }
        for (int j = 0; j <= Nx; ++j) {
            U[n + 1](0, j) = 0;
            U[n + 1](Ny, j) = 0;
        }
    }


    return U; 
}

int main ()
{
    // int N = 100; //spatial points 
    // int M = 500; //time points
    // double c = 2.0; //wave speed 
    // std::ofstream file ("solution.csv"); 
    
    // auto solution = Upwind(N, M, c); 

    // if (file)
    // {
    //     file << solution; 
        
    //     file.close(); 
    // }


    const int Nx = 200, Ny = 200, Nt = 500;
    const double c = 1.0;
    std::ofstream file ("solution.csv"); 

    auto results = Leapfrog(Nx, Ny, Nt, c); 

    if (file)
    {
        for (int i =0; i < results.size(); ++i)
        {
            for (int r = 0; r <= Ny; ++r)
            {
                for (int c = 0; c <= Nx; ++c) 
                {
                    file << results[i](r, c);
                    if (c < Nx) file << ",";  
                }
                file << "\n";  
            }
            file << "\n";
        }

        file.close(); 
    }

}