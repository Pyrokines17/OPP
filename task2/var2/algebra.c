#include <stdlib.h>
#include <math.h>
#include "algebra.h"

double scalarMul(const double* first, const double* second, int size) {
    double result = 0;
    for (int i = 0; i < size; ++i) {
        result += first[i] * second[i];
    }
    return result;
}

double* matrixMul(const double* matrix, const double* vector, int size) {
    double* result = calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    return result;
}

double* numberMul(const double* vector, double number, int size) {
    double* result = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i) {
        result[i] = number * vector[i];
    }
    return result;
}

double* sub(const double* first, const double* second, int size) {
    double* result = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i) {
        result[i] = first[i] - second[i];
    }
    return result;
}

double* sum(const double* first, const double* second, int size) {
    double* result = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i) {
        result[i] = first[i] + second[i];
    }
    return result;
}

double norm(double* vector, int size) {
    return sqrt(scalarMul(vector, vector, size));
}
