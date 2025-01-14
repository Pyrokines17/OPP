#include <stdlib.h>
#include <math.h>
#include "algebra.h"

double scalarMul(const double* first, const double* second, int size) {
    double result = 0;

    #pragma omp parallel for schedule(static) reduction(+:result)
    for (int i = 0; i < size; ++i) {
        result += first[i] * second[i];
    }

    return result;
}

double* matrixMul(const double* matrix, const double* vector, int size) {
    double* result = calloc(size, sizeof(double));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < size; ++i) {
        double buf = 0;

        for (int j = 0; j < N; ++j) {
            buf += matrix[i * size + j] * vector[j];
        }
        
        result[i] = buf;
    }

    return result;
}

double* numberMul(const double* vector, double number, int size) {
    double* result = malloc(sizeof(double) * size);

    #pragma omp for
    for (int i = 0; i < size; ++i) {
        result[i] = number * vector[i];
    }

    return result;
}

double* sub(const double* first, const double* second, int size) {
    double* result = malloc(sizeof(double) * size);

    #pragma omp for
    for (int i = 0; i < size; ++i) {
        result[i] = first[i] - second[i];
    }

    return result;
}

double* sum(const double* first, const double* second, int size) {
    double* result = malloc(sizeof(double) * size);

    #pragma omp for
    for (int i = 0; i < size; ++i) {
        result[i] = first[i] + second[i];
    }

    return result;
}

double norm(double* vector, int size) {
    return sqrt(scalarMul(vector, vector, size));
}
