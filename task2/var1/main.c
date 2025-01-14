#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "algebra.h"

#define E 0.00001

double randfrom(double min, double max) {
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main() {
    double begin = (double)clock() / CLOCKS_PER_SEC;

    double* A = malloc(sizeof(double) * N * N);

    srand(time(NULL));    
    double ra;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j <= i) {
                A[i * N + j] = randfrom(-1, 1);
                A[j * N + i] = A[i * N + j];
            }
        }
    }

    double* b = malloc(sizeof(double) * N);

    for (int i = 0; i < N; ++i) {
        b[i] = randfrom(-1, 1);
    }

    double* x = malloc(sizeof(double) * N);

    for (int i = 0; i < N; ++i) {
        x[i] = 1;
    }

    double* r = sub(b, matrixMul(A, x, N), N);
    double* z = malloc(sizeof(double) * N);
    memcpy(z, r, N * sizeof(double));
    double* r1, a, q, *tempVec, tempRes, newTempRes;
    int count = 0;
    int flag = 0;

    tempRes = scalarMul(r, r, N);

    double norma = norm(b, N);
    double div = sqrt(tempRes) / norma;

    while (div > E || flag < 3) {
        tempVec = matrixMul(A, z, N);

        a = tempRes / scalarMul(tempVec, z, N);
        x = sum(x, numberMul(z, a, N), N);
        r1 = sub(r, numberMul(tempVec, a, N), N);
        newTempRes = scalarMul(r1, r1, N);
        q = newTempRes / tempRes;
        z = sum(r1, numberMul(z, q, N), N);

        memcpy(r, r1, N * sizeof(double));
        tempRes = newTempRes;
        ++count;

        if (count > 50000) {
            printf("too much iteration's\n");
            exit(1);
        }

        div = sqrt(tempRes) / norma;

        if (div < E) {
                ++flag;
        } else {
                flag = 0;
        }
    }

    double end = (double)clock() / CLOCKS_PER_SEC;

    for (int i = 0; i < N; ++i) {
        printf("%f ", x[i]);
    }
    printf("\n%f -- %d\n", end - begin, count);

    return 0;
}
