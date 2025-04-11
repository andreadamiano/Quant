#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>

int main ()
{
    typedef Eigen::Matrix<double, 4,4> Matrix4x4; //compile time known size 
    Matrix4x4 m;

    // Eigen::MatrixXd m (4,4);  //runtime known size  

    m <<    7 , 3 , -1, 2 ,
            3 , 8 , 1 , -4,
            -1, 1 , 4 , -1,
            2 , -4, -1, 6 ;

    std::cout << "Matrix\n" << m << "\n"; 

    //LU decomposition
    Eigen::PartialPivLU<Matrix4x4> lu(m); 
    std::cout << "LU decomposition matrix\n" << lu.matrixLU() << "\n\n"; 

    //L (lower triangular matrix)
    Matrix4x4 l = Eigen::MatrixXd::Identity(4,4); //identity matrix 
    std::cout << "Identity Matrix \n" << l << "\n\n"; 

    l.triangularView<Eigen::StrictlyLower>() = lu.matrixLU(); 
    std::cout << "L matrix\n" << l << "\n\n"; 
    
    
    //U (upper triangular view)
    Matrix4x4 u = lu.matrixLU().triangularView<Eigen::Upper>(); 
    std::cout << "U matrix\n" << u << "\n\n"; 

    

}