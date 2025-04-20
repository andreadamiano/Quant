#include "Matrix.h"
#include <iostream>
#include <Eigen/Dense>


// int main ()
// {
//     Matrix matrix1(10, 10, 1.0); 
//     Matrix matrix2 (10, 10, 2.0); 

//     Matrix matrix3 = matrix1 + matrix2; 


//     for (int i =0; i< matrix3.nRows(); ++i)
//     {
//         for (int j=0; j< matrix3.nCols(); ++j)
//             std::cout << matrix3(i, j) <<" "; 
//         std::cout << "\n";
//     }
// }

int main ()
{
    Eigen::MatrixXd m(2,3); //X says that the size is dynamic and d stands for double 
    Eigen::MatrixXd m2(3,2);
    m << 1,2,3,4,5,6; 
    m2 << 10,11,12,13,14,15;
    
    auto result = m * m2; 

    std::cout << result << "\n";
    
    Eigen::Vector3d v1 (1.2, 2.2, 3.2); 
    Eigen::Vector3d v2 (4,5,6); 

    std::cout << v1 +v2 << "\n"; 

    //matrix transposition 
    Eigen::MatrixXd m3 (3,3); 
    m3 << 1,2,3,4,5,6,7,8,9; 

    std::cout << "m3\n" << m3 << "\n"; 
    std::cout << "transposed matrix:\n " << m3.transpose(); 
    std::cout << "\nm3\n" << m3 << "\n"; 

    //vector opearations 
    Eigen::Vector3d v3 (1,2,3); 
    Eigen::Vector3d v4 (4,5,6); 

    std::cout << "dot product " << v3.dot(v4) << "\n"; 
    std::cout << "cross product\n " << v3.cross(v4) << "\n"; 

    //reduction
    Eigen::MatrixXd m4(3,3); 
    m4 << 1,2,3,4,5,6,7,8,9;
    std::cout << "\nm4\n" << m4 << "\n"; 

    std::cout << "sum " << m4.sum() << "\n"; //sum of all eleemntes in the matrix 
    std::cout << "product " << m4.prod() << "\n"; //prod of all eleemntes in the matrix 
    std::cout << "mean " << m4.sum() << "\n"; //mean average of all eleemntes in the matrix 
    std::cout << "min " << m4.minCoeff() << "\n"; //min of all eleemntes in the matrix 
    std::cout << "max " << m4.maxCoeff() << "\n"; //sum of all eleemntes in the matrix 
    std::cout << "trace " << m4.trace() << "\n"; //trace of a matrix (sum of diagonal elements)




    
}