#pragma once

#define N 2016

double scalarMul(const double* first, const double* second, int size);
double* matrixMul(const double* matrix, const double* vector, int size);
double* numberMul(const double* vector, double number, int size);
double* sub(const double* first, const double* second, int size);
double* sum(const double* first, const double* second, int size);
double norm(double* vector, int size);
