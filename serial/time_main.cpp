#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <time.h>

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
    size_t i, j;

    for(i=0;i<source.size();i++)
    {
        for(j=0;j<source[i].size();j++)
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
    clock_t step_start, step_end;
    clock_t total_start = clock();
    size_t n = A.size();

    if (n == 2)
    {
        // Not printing for 2x2 to prevent millions of console prints which will freeze your terminal.
        // If you specifically want to see 2x2 times, you can add print statements here.
        return Inverse2x2(A);
    }

    // Split
    step_start = clock();
    Matrix A11, A12, A21, A22;

    SplitMatrix(A,
                A11,
                A12,
                A21,
                A22);
    step_end = clock();
    cout << "[Size " << n << "] SplitMatrix took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

    // FIRST TRY:
    // Use A11 as Pivot
    try
    {
        step_start = clock();
        Matrix invA11 = RecursiveInverse(A11);
        step_end = clock();
        cout << "[Size " << n << "] RecursiveInverse(A11) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        // Schur Complement
        step_start = clock();
        Matrix temp1 = Multiply(A21, invA11);
        Matrix temp2 = Multiply(temp1, A12);

        Matrix S = Subtract(A22, temp2);
        step_end = clock();
        cout << "[Size " << n << "] Schur Complement setup took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        step_start = clock();
        Matrix invS = RecursiveInverse(S);
        step_end = clock();
        cout << "[Size " << n << "] RecursiveInverse(S) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        // P, Q, R, T
        step_start = clock();
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
        step_end = clock();
        cout << "[Size " << n << "] P, Q, R, T block calculation took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        step_start = clock();
        Matrix result = CombineMatrix(P, Q, R, T);
        step_end = clock();
        cout << "[Size " << n << "] CombineMatrix took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";
        
        clock_t total_end = clock();
        cout << ">>> [Size " << n << "] TOTAL level execution took: " << double(total_end - total_start) / CLOCKS_PER_SEC << " s <<<\n\n";

        return result;
    }
    catch (...)
    {
        // SECOND TRY:
        // Use A22 as Pivot
        cout << "[Size " << n << "] First try failed, attempting A22 pivot...\n";
        
        step_start = clock();
        Matrix invA22 = RecursiveInverse(A22);
        step_end = clock();
        cout << "[Size " << n << "] RecursiveInverse(A22) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        // Schur Complement
        step_start = clock();
        Matrix temp1 = Multiply(A12, invA22);
        Matrix temp2 = Multiply(temp1, A21);

        Matrix S = Subtract(A11, temp2);
        step_end = clock();
        cout << "[Size " << n << "] Schur Complement setup (Try 2) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        step_start = clock();
        Matrix invS = RecursiveInverse(S);
        step_end = clock();
        cout << "[Size " << n << "] RecursiveInverse(S) (Try 2) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        step_start = clock();
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
        step_end = clock();
        cout << "[Size " << n << "] P, Q, R, T block calculation (Try 2) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";

        step_start = clock();
        Matrix result = CombineMatrix(P, Q, R, T);
        step_end = clock();
        cout << "[Size " << n << "] CombineMatrix (Try 2) took: " << double(step_end - step_start) / CLOCKS_PER_SEC << " s\n";
        
        clock_t total_end = clock();
        cout << ">>> [Size " << n << "] TOTAL level execution took: " << double(total_end - total_start) / CLOCKS_PER_SEC << " s <<<\n\n";

        return result;
    }

    throw runtime_error("Matrix is singular.");
}

int main()
{
    clock_t start_time, end_time;

   /* Matrix M =
    {
       {8, 1, 2, 3, 4, 5, 6, 7},
        {1, 8, 0, 0, 0, 0, 0, 0},
        {2, 0, 8, 0, 0, 0, 0, 0},
        {3, 0, 0, 8, 0, 0, 0, 0},
        {4, 0, 0, 0, 8, 0, 0, 0},
        {5, 0, 0, 0, 0, 8, 0, 0},
        {6, 0, 0, 0, 0, 0, 8, 0},
        {7, 0, 0, 0, 0, 0, 0, 8}
    }; */


    Matrix M = getDiagonallyDominant(1024);

    cout << "Recursive start..." << endl;

    start_time = clock();
    Matrix inv = RecursiveInverse(M);
    end_time = clock();
    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;

    cout << "Time taken to compute the inverse: " << elapsed_time << " seconds" << endl;

    //cout<<"Inverse Matrix: "<<endl;
    //PrintMatrix(inv);


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
