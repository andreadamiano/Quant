#include <Eigen/Dense>
#include <iostream>

using namespace Eigen; 

int main()
{
    //sample housing data: [size (sqft), bedrooms, price ($)]
    MatrixXd data (5, 3); 

    data << 2104, 3, 399900,
            1600, 3, 329900,
            2400, 3, 369000,
            1416, 2, 232000,
            3000, 4, 539900;

    //X, y
    MatrixXd X = data.leftCols(2);
    MatrixXd y = data.rightCols(1); 

    //scale 
    X.col(0) = X.col(0)/1000; 

    std::cout << "X: \n" << X << "\n\n"; 
    std::cout << "y: \n" << y << "\n\n"; 
    
    //add intercept 
    MatrixXd X_with_intercept (X.rows(), X.cols() + 1); 
    X_with_intercept << VectorXd::Ones(X.rows()), X; 

    std::cout << "X: \n" << X_with_intercept << "\n\n"; 

    //QR decomposition
    Eigen::HouseholderQR<MatrixXd> qr(X_with_intercept); 
    // MatrixXd Q = qr.householderQ(); 
    // MatrixXd R = qr.matrixQR().triangularView<Upper>(); 

    //predict theta
    VectorXd theta = qr.solve(y); 

    std::cout << "theta: \n" << theta << "\n\n"; 

    //make predictions 
    VectorXd preds = X_with_intercept*theta; 
    std::cout << "preds: \n" << preds << "\n\n"; 

    //compute r^2
    double y_mean = y.mean();
    double ss_tot = (y.array() - y_mean).square().sum();
    double ss_res = (y - preds).array().square().sum();
    double r_squared = 1 - (ss_res / ss_tot);
    
    std::cout << "R-squared: " << r_squared << "\n";

}

