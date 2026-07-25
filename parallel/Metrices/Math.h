#ifndef MATH_H
#define MATH_H

#include <vector>

using namespace std;

typedef vector<vector<double>> Matrix;

Matrix Add(const Matrix& A, const Matrix& B);
Matrix Subtract(const Matrix& A, const Matrix& B);
Matrix Multiply(const Matrix& A, const Matrix& B);
Matrix ScalarMultiply(const Matrix& A, double scalar);

#endif
