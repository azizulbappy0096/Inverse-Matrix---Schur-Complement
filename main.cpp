#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

#include "Metrices/Metrices.h"
#include "Metrices/Validate.h"
#include "Metrices/Math.h"

using namespace std;

typedef vector<vector<double>> Matrix;


void PrintMatrix(const Matrix& M)
{
    for(size_t i = 0; i < M.size(); i++)
    {
        for(size_t j = 0; j < M[i].size(); j++)
        {
            cout << fixed << setprecision(3) << M[i][j] << "\t";
        }

        cout << endl;
    }

    cout << endl;
}

Matrix GetBlock(const Matrix& M,
                int row,
                int col,
                int size)
{
    Matrix block(size, vector<double>(size));

    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            block[i][j]=M[row+i][col+j];
        }
    }

    return block;
}

void SetBlock(Matrix& destination,
              const Matrix& source,
              int row,
              int col)
{
    for(size_t i=0;i<source.size();i++)
    {
        for(size_t j=0;j<source[i].size();j++)
        {
            destination[row+i][col+j]=source[i][j];
        }
    }
}


void SplitMatrix(
        const Matrix& A,
        Matrix& A11,
        Matrix& A12,
        Matrix& A21,
        Matrix& A22)
{
    int half = A.size() / 2;

    A11 = GetBlock(A, 0,    0,    half);
    A12 = GetBlock(A, 0,    half, half);
    A21 = GetBlock(A, half, 0,    half);
    A22 = GetBlock(A, half, half, half);
}


Matrix CombineMatrix(
        const Matrix& P,
        const Matrix& Q,
        const Matrix& R,
        const Matrix& T)
{
    int half = P.size();

    Matrix result(
        half * 2,
        vector<double>(half * 2)
    );

    SetBlock(result, P, 0,    0);
    SetBlock(result, Q, 0,    half);
    SetBlock(result, R, half, 0);
    SetBlock(result, T, half, half);

    return result;
}



Matrix Inverse2x2(const Matrix& A)
{
    double a = A[0][0];
    double b = A[0][1];
    double c = A[1][0];
    double d = A[1][1];

    double det = (a * d) - (b * c);

    const double EPS = 1e-12;

    if (fabs(det) < EPS)
    {
        PrintMatrix(A);

        throw runtime_error("Matrix is singular (from Inverse2x2 func)");
    }

    Matrix inverse(2, vector<double>(2));

    inverse[0][0] =  d / det;
    inverse[0][1] = -b / det;
    inverse[1][0] = -c / det;
    inverse[1][1] =  a / det;

    return inverse;
}


Matrix RecursiveInverse(const Matrix& A)
{

    if (A.size() == 2)
    {
        return Inverse2x2(A);
    }


    // Split
    Matrix A11, A12, A21, A22;

    SplitMatrix(A,
                A11,
                A12,
                A21,
                A22);

    // FIRST TRY:
    // Use A11 as Pivot
    try
    {
        Matrix invA11 = RecursiveInverse(A11);

        // Schur Complement
        Matrix temp1 = Multiply(A21, invA11);
        Matrix temp2 = Multiply(temp1, A12);

        Matrix S = Subtract(A22, temp2);

        Matrix invS = RecursiveInverse(S);


        // P
        Matrix P = Multiply(invA11, A12);
        P = Multiply(P, invS);
        P = Multiply(P, A21);
        P = Multiply(P, invA11);
        P = Add(invA11, P);


        // Q
        Matrix Q = Multiply(invA11, A12);
        Q = Multiply(Q, invS);
        Q = ScalarMultiply(Q, -1.0);

        // R
        Matrix R = Multiply(invS, A21);
        R = Multiply(R, invA11);
        R = ScalarMultiply(R, -1.0);

        // T
        Matrix T = invS;

        return CombineMatrix(P, Q, R, T);
    }
    catch (...)
    {
        // SECOND TRY:
        // Use A22 as Pivot
        Matrix invA22 = RecursiveInverse(A22);


        // Schur Complement
        Matrix temp1 = Multiply(A12, invA22);
        Matrix temp2 = Multiply(temp1, A21);

        Matrix S = Subtract(A11, temp2);

        Matrix invS = RecursiveInverse(S);

        // P
        Matrix P = invS;


        // Q
        Matrix Q = Multiply(invS, A12);
        Q = Multiply(Q, invA22);
        Q = ScalarMultiply(Q, -1.0);


        // R
        Matrix R = Multiply(invA22, A21);
        R = Multiply(R, invS);
        R = ScalarMultiply(R, -1.0);


        // T
        Matrix T = Multiply(invA22, A21);
        T = Multiply(T, invS);
        T = Multiply(T, A12);
        T = Multiply(T, invA22);
        T = Add(invA22, T);

        return CombineMatrix(P, Q, R, T);
    }

    throw runtime_error("Matrix is singular.");
}


int main()
{

    /*Matrix M =
    {
       {5, 7, 9, 0},
        {4, 3, 8, 0},
        {7, 5, 6, 0},
        {0, 0, 0, 0}
    };*/

    Matrix M = getDiagonallyDominant(16);

    Matrix inv = RecursiveInverse(M);

    cout<<"Inverse Matrix: "<<endl;
    PrintMatrix(inv);


    cout<<"\n-------------------------------------------"<<endl;
    cout << "Verification: ";
    if (ValidateInverse(M, inv))
    {
        cout<<"(A * A^-1) = I"<<endl;
    }else
    {
        cout<<"(A * A^-1) != I"<<endl;
    }

    return 0;
}
