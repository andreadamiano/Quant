#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace Eigen; 

//Lapce equation: uxx + uyy = 0

// boundary conditions 
//top , right, bottom = 0 / left = 100*sin(pi*y)

constexpr double M_PI = std::acos(-1.0);

MatrixXd Jacobi (const int& N, const int& M, const double tol = 1e-6, const int max_iter = 1000)
{  
    MatrixXd U = MatrixXd::Zero(N+1, M+1); 
    MatrixXd U_new = MatrixXd::Zero(N+1, M+1); 
    std::vector<double> max_diff_history; //store convergence history
    int iter =0; 
    

    //boundary conditions
    for (int j=0; j< N+1; ++j)
        U(0, j) = U_new(0, j) = 100.0 * sin(M_PI * static_cast<double>(j)/M); 

    while(iter <max_iter) 
    {
        double max_diff = 0.0; 
        for (int i =1; i < N; ++i)
        {
            //update interior points 
            for (int j=1; j< M; ++j)
            {
                U_new(i, j) = 0.25 * (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
                max_diff = std::max(max_diff, std::abs(U_new(i,j) - U(i,j))); 
            }
        }

        max_diff_history.push_back(max_diff); 
        iter++; 

        if (max_diff < tol)
            break; 

        U= U_new;
    }

    //convergence summary
    std::cout << "Jacobi converged after " << iter << " iterations.\n";
    std::cout << "Final max error: " << max_diff_history.back() << "\n";


    return U; 
}

//update values in place 
MatrixXd GaussSeidal (const int& N, const int& M, const double tol = 1e-6, const int max_iter = 1000)
{  
    MatrixXd U = MatrixXd::Zero(N+1, M+1); 
    std::vector<double> max_diff_history; //store convergence history
    int iter =0; 
    

    //boundary conditions
    for (int j=0; j< N+1; ++j)
        U(0, j) = 100.0 * sin(M_PI * static_cast<double>(j)/M); 

    while(iter <max_iter) 
    {
        double max_diff = 0.0; 
        for (int i =1; i < N; ++i)
        {
            //update interior points 
            for (int j=1; j< M; ++j)
            {
                double old_val = U(i, j); 
                U(i, j) = 0.25 * (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
                max_diff = std::max(max_diff, std::abs(U(i,j) - old_val)); 
            }
        }

        max_diff_history.push_back(max_diff); 
        iter++; 

        if (max_diff < tol)
            break; 
    }

    //convergence summary
    std::cout << "Gauss Seidal converged after " << iter << " iterations.\n";
    std::cout << "Final max error: " << max_diff_history.back() << "\n";


    return U; 
}

//relaxation method
MatrixXd SOR (const int& N, const int& M, const double tol = 1e-6, const int max_iter = 1000)
{  
    MatrixXd U = MatrixXd::Zero(N+1, M+1); 
    std::vector<double> max_diff_history; //store convergence history
    int iter =0; 
    double h = 1.0/N; //if the grid is a square hx = hy , if is a rectangle we chose the biggest h 
    double omega = 2/(1+sin(M_PI *h));     
    

    //boundary conditions
    for (int j=0; j< N+1; ++j)
        U(0, j) = 100.0 * sin(M_PI * static_cast<double>(j)/M); 

    while(iter <max_iter) 
    {
        double max_diff = 0.0; 
        for (int i =1; i < N; ++i)
        {
            //update interior points 
            for (int j=1; j< M; ++j)
            {
                double old_val = U(i, j); 
                double gs_update = 0.25 * (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
                
                U(i, j) =  old_val * (1.0 -omega) + omega* gs_update; //relaxation
                max_diff = std::max(max_diff, std::abs(U(i,j) - old_val));
            }
        }

        max_diff_history.push_back(max_diff); 
        iter++; 

        if (max_diff < tol)
            break; 
    }

    //convergence summary
    std::cout << "SOR converged after " << iter << " iterations.\n";
    std::cout << "Final max error: " << max_diff_history.back() << "\n";


    return U; 
}

//solve linear system with LU decomposition
MatrixXd LUDecompostion (const int& N, const int& M)
{  
    MatrixXd U = MatrixXd::Zero(N+1, M+1); //solution matrix 

    int total_points = (N-1) * (M-1); 
    SparseMatrix<double> A (total_points, total_points); 
    A.reserve(VectorXi::Constant(total_points, 5)); //reserve accept an int ( i create a vector of all 5 , to reserve 5 non empty spaces for each row)
    VectorXd b (total_points); 

    //fill A and b 
    // -4*u_i,j + u_i+1,j + u_i-1,j + u_i,j+1 + u_i,j-1 = 0  , this is the linear equation for each interior point 
    for (int i =1; i < N; ++i)
        for (int j =1; j < M; ++j)
        {
            int k = (i-1)*(M-1)+(j-1); //map a 2d grid coordinates to 1d (excluding boundaries)

            //every row k contains a linear equation 
            A.insert(k, k) = -4.0; //diagonal elements 
            
            //neighbors 
            if (i > 1) A.insert(k, k - (M - 1)) = 1.0;  // u_{i-1,j}
            if (i < N - 1) A.insert(k, k + (M - 1)) = 1.0;  // u_{i+1,j}
            if (j > 1) A.insert(k, k - 1) = 1.0;  // u_{i,j-1}
            if (j < M - 1) A.insert(k, k + 1) = 1.0;  // u_{i,j+1}

            //fill b vector (which collect all the known boundary terms of the linear system)
            if (i==1)
                b(k) = -100.0 * sin(M_PI * static_cast<double>(j)/M); //the term changes since is moved to the right hand side of the equation

            else
                b(k) = 0.0; 

        }

    //solve linear system 
    SparseLU<SparseMatrix<double> > solver; 
    solver.compute(A); //LU decomposition
    VectorXd u_interior = solver.solve(b); 

    //reconstruct full solution
    for (int j=0; j< M+1; ++j)
        U(0, j) = 100.0 * sin(M_PI * static_cast<double>(j)/M); 

    for (int i = 1; i < N; ++i) 
        for (int j = 1; j < M; ++j)
        {
            int k = (i-1)* (M-1) + (j-1); 
            U(i, j) = u_interior(k); //fill interior points 
        }
    
    return U; 
}


//Lapce equation: uxx + uyy = 0

// boundary conditions:
//left = 100*sin(pi*y)

//Neumann boundary conditions
//right: ux(1, y) =0
//top: uy(x,1) =0
//bottom: uy(x,0) =0

MatrixXd SOR_2 (const int& N, const int& M, const double tol = 1e-6, const int max_iter = 1000)
{  
    MatrixXd U = MatrixXd::Zero(N+1, M+1); 
    std::vector<double> max_diff_history; //store convergence history
    int iter =0; 
    double h = 1.0/N; //if the grid is a square hx = hy , if is a rectangle we chose the biggest h 
    double omega = 2/(1+sin(M_PI *h));     
    

    //left boundary conditions (Dirichelet)
    // for (int i=0; i< N+1; ++i)
    // {
    //     U(i, 0) = 1000.0 * sin(M_PI * static_cast<double>(i)/N); 
    //     U(i, M) = 1000.0 * sin(M_PI * static_cast<double>(i)/N); 
    // }

    // for (int j=0; j< N+1; ++j)
    // {
    //     U(0, j) = 1000.0 * sin(M_PI * static_cast<double>(j)/N); 
    //     U(N, j) = 1000.0 * sin(M_PI * static_cast<double>(j)/N); 
    // }

    for (int i=0; i< N+1; ++i)
        U(i, 0) = 1000.0 * sin(M_PI * static_cast<double>(i)/N); 

    while(iter <max_iter) 
    {
        double max_diff = 0.0; 
        for (int i =1; i < N; ++i)
        {
            //update interior points 
            for (int j=1; j< M; ++j)
            {
                double old_val = U(i, j); 
                double gs_update = 0.25 * (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
                
                U(i, j) =  old_val * (1.0 -omega) + omega* gs_update; //relaxation
                max_diff = std::max(max_diff, std::abs(U(i,j) - old_val));
            }
        }

        //Neumann boundary conditions 

        //one sided approximation
        // for (int j = 0; j <= M; ++j)
        //     U(N, j) = U(N-1, j);  // right (ux=0)

        // for (int i = 0; i <= N; ++i)
        //     U(i, M) = U(i, M-1);  // top (uy=0)

        // for (int i = 1; i <= N; ++i)
        //     U(i, 0) = U(i, 1);    // bottom (uy=0), skip i=0

        //ghost point approximation
        //right boundary 
        for (int j = 1; j < M; ++j) {
            double old_val = U(N, j);
            double gs_update = (2.0 * U(N-1, j) + U(N, j+1) + U(N, j-1)) / 4.0;
            U(N, j) = old_val * (1.0 - omega) + omega * gs_update;
            max_diff = std::max(max_diff, std::abs(U(N, j) - old_val));
        }

        //bottom boundary
        for (int i = 1; i < N; ++i) {
            double old_val = U(i, 0);
            double gs_update = (U(i+1, 0) + U(i-1, 0) + 2.0 * U(i, 1)) / 4.0;
            U(i, 0) = old_val * (1.0 - omega) + omega * gs_update;
            max_diff = std::max(max_diff, std::abs(U(i, 0) - old_val));
        }

        //top boundary 
        for (int i = 1; i < N; ++i) {
            double old_val = U(i, M);
            double gs_update = (U(i+1, M) + U(i-1, M) + 2.0 * U(i, M-1)) / 4.0;
            U(i, M) = old_val * (1.0 - omega) + omega * gs_update;
            max_diff = std::max(max_diff, std::abs(U(i, M) - old_val));
        }

        //right-bottom corner (righ + bottom Neumann conditions)
        double old_val_corner = U(N, 0);
        double gs_corner = (U(N-1, 0) + U(N, 1)) / 2.0;
        U(N, 0) = old_val_corner * (1.0 - omega) + omega * gs_corner;
        max_diff = std::max(max_diff, std::abs(U(N, 0) - old_val_corner));

        //right-top corner (righ + top Neumann conditions)
        old_val_corner = U(N, M);
        gs_corner = (U(N-1, M) + U(N, M-1)) / 2.0;
        U(N, M) = old_val_corner * (1.0 - omega) + omega * gs_corner;
        max_diff = std::max(max_diff, std::abs(U(N, M) - old_val_corner));

        
        max_diff_history.push_back(max_diff); 
        iter++; 

        if (max_diff < tol)
            break; 
    }

    //convergence summary
    std::cout << "SOR converged after " << iter << " iterations.\n";
    std::cout << "Final max error: " << max_diff_history.back() << "\n";


    return U; 
}

int main ()
{
    int N = 200; 
    int M = 200; 
    std::ofstream file ("solution.csv"); 

    MatrixXd solution = SOR_2(N, M); 

    if (file)
    {
        file << solution; 
        file.close(); 
    }
}