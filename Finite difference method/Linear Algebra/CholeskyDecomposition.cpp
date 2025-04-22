#include <Eigen/Dense>
#include <iostream>

using namespace Eigen; 

int main ()
{
    typedef Matrix<double, 4, 4> Matrix4x4; //compile time matrix dimension

    Matrix4x4 m; 
    m <<  6 , 3 , 4 , 8 ,
        3 , 6 , 5 , 1 ,
        4 , 5 , 10 , 7 ,
        8 , 1 , 7 , 25;

    Eigen:: LLT<Matrix4x4> llt(m); 

    Matrix4x4 l = llt.matrixL(); 

    std::cout << "L: \n" << l << "\n\n"; 

    //check L*L^T
    std::cout << "m: \n" << l * l.transpose(); 

}