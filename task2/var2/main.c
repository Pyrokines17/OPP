#include <stdlib.h>
#include <mpi.h>
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

int main(int argc, char* argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    double begin = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *A, *b = malloc(sizeof(double) * N),
        *x = malloc(sizeof(double) * N);
    int part = N / size;
    double* A1 = malloc(sizeof(double) * N * part);

    if (rank == 0) {
        A = malloc(sizeof(double) * N * N);
        srand(time(NULL));

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (j <= i) {
                A[i * N + j] = randfrom(-1, 1);
                A[j * N + i] = A[i * N + j];
            }
            }
        }

        for (int i = 0; i < N; ++i) {
            b[i] = randfrom(-1, 1);
        }

        for (int i = 0; i < N; ++i) {
            x[i] = 1;
        }
    }

    MPI_Scatter(A, N * part, MPI_DOUBLE, A1, N * part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

double* partOfMul = matrixMul(A1, x, part),
        *r;

    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            r = sub(b + i * part, partOfMul, part);
        }
    }

    double* z = malloc(sizeof(double) * part);
    memcpy(z, r, sizeof(double) * part);
    double a, q, tempRes, tempRes1, tempRes2, *tempVec, up, up1, down;
    double* zFull = malloc(sizeof(double) * N),
        *resPart = malloc(sizeof(double) * part),
        *r1;
    int count = 0,
        flag = 0;

    MPI_Scatter(x, part, MPI_DOUBLE, resPart, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    tempRes = scalarMul(r, r, part);
    MPI_Allreduce(&tempRes, &up, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double norma = norm(b, N);
    double div = sqrt(up) / norma;

while (div > E || flag < 3) {
        MPI_Allgather(z, part, MPI_DOUBLE, zFull, part, MPI_DOUBLE, MPI_COMM_WORLD);

        tempVec = matrixMul(A1, zFull, part);
        tempRes1 = scalarMul(tempVec, z, part);

        MPI_Allreduce(&tempRes1, &down, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        a = up / down;
        resPart = sum(resPart, numberMul(z, a, part), part);
        r1 = sub(r, numberMul(tempVec, a, part), part);
        tempRes2 = scalarMul(r1, r1, part);

        MPI_Allreduce(&tempRes2, &up1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        q = up1 / up;
        z = sum(r1, numberMul(z, q, part), part);

        memcpy(r, r1, sizeof(double) * part);
        up = up1;
        ++count;

        if (count > 50000) {
            if (rank == 0) {
                printf("too much iteration's\n");
            }
            exit(1);
        }

        div = sqrt(up) / norma;

        if (div < E) {
                ++flag;
        } else {
                flag = 0;
        }
    }

MPI_Gather(resPart, part, MPI_DOUBLE, x, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double end = MPI_Wtime();

    if (rank == 0) {
        for (int i = 0; i < N; ++i) {
            printf("%f ", x[i]);
        }
        printf("\n%f -- %d\n", end - begin, count);
    }

    MPI_Finalize();
    return 0;
}