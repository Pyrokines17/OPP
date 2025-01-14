#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define sizeM 48

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    double begin = MPI_Wtime();

    int p1 = atoi(argv[1]), p2 = atoi(argv[2]);

    int myId;
    int n1 = sizeM, 
        n2 = sizeM, 
        n3 = sizeM; 

    MPI_Comm gorComm, vertComm;

    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Comm_split(MPI_COMM_WORLD, myId/p2, myId, &gorComm);
    MPI_Comm_split(MPI_COMM_WORLD, myId%p2, myId, &vertComm);

    int partA = n1*n2/p1,
        partB = n2*n3/p2;

    double *A, *B, *res,
            *A1 = (double*)malloc(partA*sizeof(double)), 
            *B1 = (double*)malloc(partB*sizeof(double));

    if (myId == 0) {
        A = (double*)malloc(n1*n2*sizeof(double));
        B = (double*)malloc(n2*n3*sizeof(double));

        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                A[i * n2 + j] = i * n2 + j;
            }
        }

        double maxV = n2 * n3 - 1;

        for (int i = 0; i < n2; ++i) {
            for (int j = 0; j < n3; ++j) {
                B[i * n3 + j] = maxV - (i*n3+j);
            }
        }
    }

    if (myId%p2 == 0) {
        MPI_Scatter(A, partA, MPI_DOUBLE, A1, partA, MPI_DOUBLE, 0, vertComm);
    }

    MPI_Bcast(A1, partA, MPI_DOUBLE, 0, gorComm);

    MPI_Datatype col, colType;

    MPI_Type_vector( n2 , n3/p2 , n3 , MPI_DOUBLE , &col);
    MPI_Type_commit( &col);

    MPI_Type_create_resized(col, 0, sizeof(double) * partB/n2, &colType);
    MPI_Type_commit(&colType);

    if (myId/p2 == 0) {
        MPI_Scatter(B, 1, colType, B1, partB, MPI_DOUBLE, 0, gorComm);
    }

    MPI_Bcast(B1, partB, MPI_DOUBLE, 0, vertComm);

    int weiLocRes = partB/n2, heiLocRes = partA/n2;
    double* locRes = calloc(weiLocRes*heiLocRes, sizeof(double));

    for (int k = 0; k < heiLocRes; ++k) {
        for (int i = 0; i < n2; ++i) {
            for (int j = 0; j < weiLocRes; ++j) {
                locRes[k * weiLocRes + j] += A1[k*n2+i] * B1[i*weiLocRes+j];
            }
        }
    }

    if (myId == 0) {
        res = malloc(sizeof(double) * n1 * n3);

        for (int i = 0; i < heiLocRes; ++i) {
            for (int j = 0; j < weiLocRes; ++j) {
                res[i * n3 + j] = locRes[i * weiLocRes + j];
            }
        }
    } else {
        MPI_Send(locRes, heiLocRes*weiLocRes, MPI_DOUBLE, 0, myId, MPI_COMM_WORLD);
    }

    int shift;
    MPI_Datatype block;
    MPI_Type_vector(heiLocRes, weiLocRes, n3, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);

    if (myId == 0) {
        for (int i = 1; i < p1*p2; ++i) {
            shift = (i/p2) * partA + (i%p2) * weiLocRes;
            MPI_Recv(res + shift, 1, block, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        free(A); free(B); free(res);
    }

    free(A1); free(B1);
    MPI_Type_free(&col);
    MPI_Type_free(&block);

    double end = MPI_Wtime();

    if (myId == 0) {
        printf("Time: %f\n", end-begin);
    }

    MPI_Finalize();
    return 0;
}
