#include <Eigen/Dense>
#include <iostream>

using namespace Eigen; 

int main()
{
    int n = 5; 
    typedef Vector<double, 5> Vector5d;
    typedef Matrix<double, 5, 5> Matrix5x5; 

    //define diagonal vectors 
    Vector5d mainDiagonal = Vector5d::Constant(n, 2); 
    Vector4d subDiagonal = Vector4d::Constant(n-1, -1); 
    Vector4d superDiagonal = Vector4d::Constant(n-1, -1); 

    //create tridiagonal matrix 
    Matrix5x5 tridiagonal = mainDiagonal.asDiagonal(); 
    tridiagonal.diagonal(-1) = subDiagonal; 
    tridiagonal.diagonal(1) = superDiagonal; 

    //modificate the element that need to change 
    tridiagonal(0,0) = 1; 
    tridiagonal(n-1,n-1) = 1; 

    tridiagonal(n-1, n-2) = 0; //last elemement of lower diagonal 
    tridiagonal(0,1) =0; //first element of upper diagonal 



    std::cout << "tridiagonal matrix\n" << tridiagonal; 
    
}