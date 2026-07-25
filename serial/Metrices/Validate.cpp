#include <cmath>
#include "Validate.h"
#include "Math.h"
#include "Metrices.h"

bool ValidateInverse(const Matrix& M, const Matrix& Inv) {
    double tolerance = 1e-9;

    Matrix MxInv = Multiply(M, Inv);
    Matrix I = getIdentity(M.size());

    for (size_t i = 0; i < MxInv.size(); i++) {
        for (size_t j = 0; j < MxInv[i].size(); j++) {
            if (fabs(MxInv[i][j] - I[i][j]) > tolerance) {
                return false;
            }
        }
    }

    return true;
}
