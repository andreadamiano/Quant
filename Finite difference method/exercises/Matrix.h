#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

template <typename T = double>
class Matrix 
{
    private:
        std::vector<std::vector <T> > matrix; 
        size_t rows; 
        size_t cols; 

    public:
        //constructor destructor 
        Matrix(size_t _rows, size_t _cols, T _initial); //parameter constructor 
        Matrix(const Matrix<T>& other ); //copy constructor 
        virtual ~Matrix() {}

        //asssignment operator 
        Matrix<T>& operator=(const Matrix<T>& other ); 
        

        //matrix operations 
        Matrix<T> operator +(const Matrix<T>& other); 
        Matrix<T>& operator+= (const Matrix<T>& other); 
        Matrix<T> operator- (const Matrix<T>& other); 
        Matrix<T>& operator-= (const Matrix<T>& other); 
        Matrix<T> operator* (const Matrix<T>& other); 
        Matrix<T>& operator*= (const Matrix<T>& other); 
        Matrix<T> transpose(); 

        //scalar-matrix operations
        Matrix<T> operator+(const T& scalar); 
        Matrix<T> operator-(const T& scalar); 
        Matrix<T> operator*(const T& scalar); 
        Matrix<T> operator/(const T& scalar); 

        //matrix-vector operations
        std::vector<T> operator*(const std::vector<T>& vector); 
        std::vector<T> diagonalVector(); 

        //access elements 
        T& operator() (const size_t& row, const size_t& col); 
        const T& operator() (const size_t& row, const size_t& col) const; 
        
        //col and row sizes
        size_t nRows() const {return rows; }
        size_t nCols() const {return cols; }

}; 


template <typename T>
Matrix<T>::Matrix(size_t _rows, size_t _cols, T _initial) : rows(_rows), cols(_cols)
{
    matrix.resize(rows); 

    for (size_t i=0; i<rows; ++i)
    {
        matrix[i].resize(cols, _initial); 
    }
}

template <typename T> 
Matrix<T>::Matrix(const Matrix<T>& other ) : rows(other.rows), cols (other.rows), matrix(other.matrix) {}


template <typename T> 
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other)
{
    if(this != &other)
    {
        rows = other.rows; 
        cols = other.cols; 

        //resize old matrix 
        matrix.resize(rows); 
        for (size_t i=0; i< rows; ++i)
            matrix[i].resize(cols); 
        
        //copy elements 
        for (size_t i =0; i< rows ; ++i)
        {
            for (size_t j=0; j< cols; ++j)
                matrix[i][j] = other(i, j); 
        }
    }

    return *this; 
}

template <typename T> 
Matrix<T> Matrix<T>::operator +(const Matrix<T>& other)
{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Error, the 2 matrix have different dimensions"); 

    Matrix<T> result(rows, cols, 0.0);
    
    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] + other(i, j); 
    }

    return result; 

}

template <typename T> 
Matrix<T>& Matrix<T>::operator +=(const Matrix<T>& other)
{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Error, the 2 matrix have different dimensions"); 
    
    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            this->matrix[i][j] += other(i, j); 
    }

    return *this; 

}

template <typename T> 
Matrix<T> Matrix<T>::operator -(const Matrix<T>& other)
{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Error, the 2 matrix have different dimensions"); 
        
    Matrix<T> result(rows, cols, 0.0);
    
    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] - other(i, j); 
    }

    return result; 

}

template <typename T> 
Matrix<T>& Matrix<T>::operator -=(const Matrix<T>& other)
{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Error, the 2 matrix have different dimensions"); 
    
    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            this->matrix[i][j] -= other(i, j); 
    }

    return *this; 

}

template <typename T> 
Matrix<T> Matrix<T>::operator *(const Matrix<T>& other)
{
    if (cols != other.rows)
        throw std::runtime_error("Error, the inner dimension are not equals"); 
        
    Matrix<T> result(rows, other.cols, 0.0); //outer dimensions 
    
    for (size_t i =0; i< result.rows ; ++i) 
    {
        for (size_t k =0; k< cols; ++k) //swapped j with k to improve cache locality 
        {
            for (size_t j=0; j< result.cols; ++j)
                result(i,j) += this->matrix[i][k] * other(k,j); 
        }
    }

    return result; 

}

template <typename T> 
Matrix<T>& Matrix<T>::operator *=(const Matrix<T>& other)
{
    if (cols != other.rows)
        throw std::runtime_error("Error, the inner dimension are not equals"); 
        
    Matrix<T> result = (*this) * other; 
    *this = result; 

    return *this; 

}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T& scalar)
{
    Matrix<T> result(rows, cols, 0.0); 

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] + scalar;  
    }

    return result;  
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const T& scalar)
{
        Matrix<T> result(rows, cols, 0.0);

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] - scalar;  
    }

    return result;  
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar)
{

    Matrix<T> result(rows, cols, 0.0);

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] * scalar;  
    }

    return result;  
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const T& scalar)
{

    Matrix<T> result(rows, cols, 0.0);

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = this->matrix[i][j] / scalar;  
    }

    return result;  
}

template <typename T> 
std::vector<T> Matrix<T>::operator*(const std::vector<T>& vector)
{
    if (cols != vector.size())
        throw std::runtime_error("Error, the inner dimension are not equals"); 

    std::vector<T> result(rows, 0.0); 

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result[i] += this->matrix[i][j] * vector[j]; 
    }

    return result; 
    

}

template <typename T>
std::vector<T> Matrix<T>::diagonalVector()
{
    if (rows != cols)
        throw std::runtime_error("Error, the matrix is not cubic"); 
    std::vector<T> result(rows, 0.0); 

    for (size_t i=0; i<<rows; ++i)
        result[i] = matrix[i][i]; 

    return result; 
}


template <typename T>
T& Matrix<T>::operator() (const size_t& row, const size_t& col)
{
    return matrix[row][col]; 
}

template <typename T>
Matrix<T> Matrix<T>::transpose()
{
    Matrix<T> result (rows,cols, 0.0); 

    for (size_t i =0; i< rows ; ++i)
    {
        for (size_t j=0; j< cols; ++j)
            result(i, j) = matrix[j][i];  
    }

    return result; 
}

template <typename T>
const T& Matrix<T>::operator() (const size_t& row, const size_t& col) const
{
    return matrix[row][col]; 
}


#endif 