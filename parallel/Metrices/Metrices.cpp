#include "Metrices.h"

// 8×8 Diagonally Dominant Matrix
Matrix get8x8Invertible() {
    return {
        {8, 1, 2, 3, 4, 5, 6, 7},
        {1, 8, 0, 0, 0, 0, 0, 0},
        {2, 0, 8, 0, 0, 0, 0, 0},
        {3, 0, 0, 8, 0, 0, 0, 0},
        {4, 0, 0, 0, 8, 0, 0, 0},
        {5, 0, 0, 0, 0, 8, 0, 0},
        {6, 0, 0, 0, 0, 0, 8, 0},
        {7, 0, 0, 0, 0, 0, 0, 8}
    };
}

// 16×16 Diagonally Dominant Matrix
Matrix get16x16Invertible() {
    return {
        {16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
        {1, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {2, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {3, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {4, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {6, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {7, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0},
        {8, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0},
        {9, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0},
        {10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0},
        {11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0},
        {12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0},
        {13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0},
        {14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0},
        {15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16}
    };
}

// Hilbert Matrix (always invertible)
Matrix get4x4Hilbert() {
    Matrix H(4, vector<double>(4));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
    return H;
}

Matrix get16x16Hilbert() {
    Matrix H(16, vector<double>(16));
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
    return H;
}

// Identity Matrix of any size
Matrix getIdentity(int n) {
    Matrix I(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.0;
    }
    return I;
}

// Generate custom diagonally dominant matrix of any size
Matrix getDiagonallyDominant(int n) {
    Matrix M(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        M[i][i] = n;  // Diagonal = n
        for (int j = i + 1; j < n; j++) {
            M[i][j] = (j - i) % 5 + 1;  // Small positive values 1-5
            M[j][i] = (i + 1) % 5 + 1;  // Small positive values 1-5
        }
    }
    return M;
}

// Simple 2×2 invertible matrix (for testing)
Matrix get2x2Invertible() {
    return {
        {2.0, 1.0},
        {1.0, 3.0}
    };
}
