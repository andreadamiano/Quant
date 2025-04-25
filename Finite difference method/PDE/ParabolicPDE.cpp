#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>          // Required for sparse matrices
#include <Eigen/IterativeLinearSolvers>  // For ConjugateGradient
#include <cmath>
#include <random>
#define M_PI 3.14159265358979323846

//PDE: ut = a*uxx  x(0,1), t(0,1)

//initial conditions: u(x, 0) = sin (pi*x)
//boundary conditions: u(0,t) = u(1, t) =0


using namespace Eigen; 

MatrixXd Crank_nicolson_LU(const int& N, const int& M , const double& alpha)
{
    MatrixXd U = MatrixXd::Zero(N+1, M+1); //sapce x time 
    double dt = 8.0/M; 
    double dx = 1.0/N; 
    double r = alpha*dt/(2*dx*dx); 

    //initial conditions
    for (int i=0; i< N; ++i)
        U(i, 0) = sin(M_PI * i*dx); 

    //boundary conditions (already taken into acocount when intialized matrix to 0s)

    //fill A 
    MatrixXd A = MatrixXd::Zero(N-1, N-1); //only internal points (space * space)

    VectorXd mainDiagonal = VectorXd::Constant(N-1, 1.0 + 2.0* r); 
    VectorXd upperDiagonal = VectorXd::Constant(N-2, -r); 

    A.diagonal(0) = mainDiagonal; 
    A.diagonal(-1) = A.diagonal(1) = upperDiagonal;  

    PartialPivLU<MatrixXd> lu(A); 

    //solve each linear system with LU deomposition
    for (int j=0; j<M; ++j) //time stepping loop 
    {
        VectorXd b (N-1); 

        //fill vector b 
        for (int i=1; i<N; ++i)
            b(i-1) = r* U(i-1, j) +(1.0 - 2.0*r) * U(i, j) + r * U(i + 1 ,j);  //b(i-1) is used to match 0 index order 

        //solve linear system 
        VectorXd u_inner = lu.solve(b); 

        //update solution
        for(int i =1; i<N; i++)
            U(i, j+1) = u_inner(i-1); 
     
    }


    return U; 
}

MatrixXd CrankNicolson_Multigrid (const int& M, const int& N, const double& alpha)
{
    MatrixXd U = MatrixXd::Zero(N+1, M+1); //sapce x time 
    double dt = 8.0/M; 
    double dx = 1.0/N; 
    double r = alpha*dt/(2*dx*dx); 

    //initial conditions
    for (int i=0; i< N; ++i)
        U(i, 0) = sin(M_PI * i*dx); 

    //boundary conditions (already taken into acocount when intialized matrix to 0s)

    //fill A (sparse matrix)
    Eigen::SparseMatrix<double> A(N-1, N-1);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3*(N-1) - 2); // pre allocate memory for 3 diagonals 

    for (int i = 0; i < N-1; ++i) {
        triplets.emplace_back(i, i, 1.0 + 2.0 * r); //diagonal
        if (i > 0)
            triplets.emplace_back(i, i-1, -r);      //lower diagonal
        if (i < N-2)
            triplets.emplace_back(i, i+1, -r);      //upper diagonal
    }
    A.setFromTriplets(triplets.begin(), triplets.end());


    //multigrid solver 
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                             Eigen::Lower|Eigen::Upper, 
                             Eigen::DiagonalPreconditioner<double>> solver;
    solver.compute(A);

    //solve each linear system with Multigrid
    for (int j=0; j<M; ++j) //time stepping loop 
    {
        VectorXd b (N-1); 

        //fill vector b 
        for (int i=1; i<N; ++i)
            b(i-1) = r* U(i-1, j) +(1.0 - 2.0*r) * U(i, j) + r * U(i + 1 ,j);  //b(i-1) is used to match 0 index order 

        //solve linear system 
        Eigen::VectorXd u_inner = solver.solve(b); 

        //update solution
        for(int i =1; i<N; i++)
            U(i, j+1) = u_inner(i-1); 
     
    }

    return U; 
}



int main ()
{
    const int M = 127; // Number of space steps (for mulitgrid must be 2^n -1)
    const int N = 1000; // Number of time steps
    const double alpha = 0.01; // Diffusion coefficient
    std::ofstream file("solution.csv"); 

    // MatrixXd solution = Crank_nicolson_LU(M, N, alpha);
    MatrixXd solution = CrankNicolson_Multigrid(M, N,alpha); 

    if (file)
    {
        file << solution; 

        file.close();
    }
}