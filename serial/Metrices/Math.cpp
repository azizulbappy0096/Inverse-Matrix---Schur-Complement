#include "Math.h"


Matrix Add(const Matrix& A,
           const Matrix& B)
{
    int rows = A.size();
    int cols = A[0].size();

    Matrix C(rows, vector<double>(cols));

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return C;
}


Matrix Subtract(const Matrix& A,
                const Matrix& B)
{
    int rows = A.size();
    int cols = A[0].size();

    Matrix C(rows, vector<double>(cols));

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }

    return C;
}


Matrix Multiply(const Matrix& A,
                const Matrix& B)
{
    int rows = A.size();
    int cols = B[0].size();
    int inner = A[0].size();


    Matrix C(rows, vector<double>(cols, 0));

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            for(int k = 0; k < inner; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}


Matrix ScalarMultiply(const Matrix& A, double scalar)
{
    Matrix result = A;

    for(size_t i = 0; i < A.size(); i++)
    {
        for(size_t j = 0; j < A[i].size(); j++)
        {
            result[i][j] *= scalar;
        }
    }

    return result;
}
