#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 115920

int main(int argc, char* argv[]) {
	int* a = malloc(sizeof(int) * N);
	int* b = malloc(sizeof(int) * N);
	long double s = 0;
	double time;

	clock_t begin = clock();

	for (int i = 0; i < N; ++i) {
		a[i] = 3;
	}

	memcpy(b, a, sizeof(int) * N);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			s += a[i] * b[j];
		}
	}

	clock_t end = clock();

	time = (end - begin) / CLOCKS_PER_SEC;

	printf("%Lf -- %f", s, time);

	return 0;
}