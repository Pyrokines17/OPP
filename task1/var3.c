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

	double start, finish;
	long double res, s = 0;

	int part = N / size;

	int* a1,* a = malloc(sizeof(int) * part);
	int* b = malloc(sizeof(int) * N);

	if (rank == 0) {
		start = MPI_Wtime();

		a1 = malloc(sizeof(int) * N);

		for (int i = 0; i < N; ++i) {
			a1[i] = 3;
		}

		memcpy(b, a1, sizeof(int) * N);
	}

	MPI_Scatter(a1, part, MPI_INT, a, part, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < part; ++i) {
		for (int j = 0; j < N; ++j) {
			s += a[i] * b[j];
		}
	}

	MPI_Reduce(&s, &res, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		finish = MPI_Wtime();

		printf("%Lf -- %f", res, finish - start);
	}

	MPI_Finalize();

	return 0;
}