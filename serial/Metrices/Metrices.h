#ifndef MATRICES_H
#define MATRICES_H

#include <vector>

using namespace std;

typedef vector<vector<double>> Matrix;

// Function declarations
Matrix get8x8Invertible();
Matrix get16x16Invertible();
Matrix get4x4Hilbert();
Matrix get16x16Hilbert();
Matrix getIdentity(int n);
Matrix getDiagonallyDominant(int n);
Matrix get2x2Invertible();

#endif
