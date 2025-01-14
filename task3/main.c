#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define E 0.00001
#define N 2016

double randfrom(double min, double max) {
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main() {
    double begin = omp_get_wtime();

    double* r = malloc(sizeof(double) * N),
            * z = malloc(sizeof(double) * N);
    double* A = malloc(sizeof(double) * N * N), 
            * b = malloc(sizeof(double) * N), 
            * x = malloc(sizeof(double) * N);
    double* r1 = malloc(sizeof(double) * N), 
            * buf = calloc(N, sizeof(double)),
            *tempVec = calloc(N, sizeof(double)),
            a, q, tempRes, newTempRes;
    int count = 0, flag = 0;
    double norma, div, buf1, down;
    srand(time(NULL));

    #pragma omp parallel
    {

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (j <= i) {
                    A[i * N + j] = randfrom(-1, 1);
                    A[j * N + i] = A[i * N + j];
                }
            }
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            b[i] = randfrom(-1, 1);
            x[i] = 1;
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            double bufVar = 0;
            for (int j = 0; j < N; ++j) {
                bufVar += A[i * N + j] * x[j];
            }  
            buf[i] = bufVar;
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            r[i] = b[i] - buf[i];
            z[i] = r[i];
        }

        #pragma omp for schedule(static) reduction(+:tempRes) 
        for (int i = 0; i < N; ++i) {
            tempRes += r[i] * r[i];
        }

        #pragma omp for schedule(static) reduction(+:buf1) 
        for (int i = 0; i < N; ++i) {
            buf1 += b[i] * b[i];
        }

        #pragma omp single
        {
            down = buf1;
            div = tempRes / down;
            buf1 = 0;
        }

        while (div > E * E || flag < 3) {
            #pragma omp for schedule(static)
            for (int i = 0; i < N; ++i) {
                double bufVar = 0;
                for (int j = 0; j < N; ++j) {
                    bufVar += A[i * N + j] * z[j];
                }  
                tempVec[i] = bufVar;
            }

            #pragma omp for schedule(static) reduction(+:buf1)
            for (int i = 0; i < N; ++i) {
                buf1 += tempVec[i] * z[i];
            }

            #pragma omp single
            {
                a = tempRes / buf1;
                newTempRes = 0;
                buf1 = 0;
            }

            #pragma omp for schedule(static)
            for (int i = 0; i < N; ++i) {
                x[i] = x[i] + z[i] * a;
                r1[i] = r[i] - tempVec[i] * a;
            }

            #pragma omp for schedule(static) reduction(+:newTempRes)
            for (int i = 0; i < N; ++i) {
                newTempRes += r1[i] * r1[i];
            }
            
            #pragma omp single
            {
                q = newTempRes / tempRes;
            }

            #pragma omp for schedule(static)
            for( int i = 0; i < N; ++i) {
                z[i] = r1[i] + z[i] * q;
            }

            #pragma omp single
            {
                memcpy(r, r1, N * sizeof(double));
                tempRes = newTempRes;
                ++count;
            }

            if (count > 50000) {
                printf("too much iteration's\n");
                exit(1);
            }

            #pragma omp single
            {
                div = tempRes / down;

                if (div < E * E) {
                    ++flag;
                } else {
                    flag = 0;
                }
            }

        }
    }

    free(r); free(z); free(A);
    free(b); free(x); free(r1);

    double end = omp_get_wtime();

    printf("%f -- %d\n", end - begin, count);

    return 0;
}