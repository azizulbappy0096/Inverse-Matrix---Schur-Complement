#ifndef VALIDATE_H
#define VALIDATE_H

#include <vector>

using namespace std;

typedef vector<vector<double>> Matrix;

bool ValidateInverse(const Matrix& M, const Matrix& Inv);

#endif
