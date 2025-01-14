#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define N 115920

int main(int argc, char* argv[]) {
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		int* as = malloc(sizeof(int) * N);
		int* bs = malloc(sizeof(int) * N);
		long double ss = 0, sbuf;
		double start, finish;

		start = MPI_Wtime();

		for (int i = 0; i < N; ++i) {
			as[i] = 3;
		}

		memcpy(bs, as, sizeof(int) * N);

		int part = N / (size - 1);

		for (int i = 1; i < size; ++i) {
			MPI_Send(bs, N, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(as + ((part * i) % N), part, MPI_INT, i, 1, MPI_COMM_WORLD);
		}

		for (int i = 1; i < size; ++i) {
			MPI_Recv(&sbuf, 1, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			ss += sbuf;
		}

		finish = MPI_Wtime();

		printf("%Lf -- %f", ss, finish - start);
	}

	if (rank != 0) {
		int part = N / (size - 1);
		int* ar = malloc(sizeof(int) * part);
		int* br = malloc(sizeof(int) * N);
		long double sr = 0;

		MPI_Recv(br, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(ar, part, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for (int i = 0; i < part; ++i) {
			for (int j = 0; j < N; ++j) {
				sr += ar[i] * br[j];
			}
		}

		MPI_Send(&sr, 1, MPI_LONG_DOUBLE, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;
}