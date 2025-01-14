#include <algorithm>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <atomic>
#include <math.h>
#include <mpi.h>
#include <queue>

#define L 100
#define MIN 25
#define BORDER 5

int rank, size;
double procRes = 0;
std::atomic<bool> othHave = false, iHave = false;

std::queue<int> tasks;
pthread_t reqThread, handleThread;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

const int toHandler = 0, toRequester = 1;

void* request_tasks(void* arg) {
    bool* havent = new bool[size]{false};
    
    do {
        pthread_mutex_lock(&mutex);
        int mySize = tasks.size();
        pthread_mutex_unlock(&mutex);

        if (mySize < MIN) {
            int countHavent = 0;

            for (int i = 1; i < size; ++i) {
                int nextRank = (rank+i)%size;
                int newTask;

                if (havent[nextRank]) {
                    ++countHavent;
                    continue;
                }

                MPI_Send(&mySize, 1, MPI_INT, nextRank, toHandler, MPI_COMM_WORLD);
                printf("Rank: %d -- send size %d to %d\n", rank, mySize, nextRank);

                MPI_Recv(&newTask, 1, MPI_INT, nextRank, toRequester, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("Rank: %d -- recv newTask %d from %d\n", rank, newTask, nextRank);

                if (newTask < 0) {
                    ++countHavent;

                    if (newTask == -2) {
                        havent[nextRank] = true;
                    }
                } else {
                    pthread_mutex_lock(&mutex);
                    tasks.push(newTask);
                    pthread_mutex_unlock(&mutex);
                    iHave = true;
                }
            }

            printf("Rank: %d -- countHavent %d\n", rank, countHavent);

            if (countHavent == size-1) {
                othHave = false;
            }
        }
    } while (othHave);

    printf("Rank: %d -- end of request\n", rank);
}

void* handle_requests(void* arg) {
    int recVar, notCount = 0;
    MPI_Status status;
    bool end = false;

    do {
        MPI_Recv(&recVar, 1, MPI_INT, MPI_ANY_SOURCE, toHandler, MPI_COMM_WORLD, &status);
        printf("Rank: %d -- recv size %d from %d\n", rank, recVar, status.MPI_SOURCE);
        int senderRank = status.MPI_SOURCE;

        if (end) {
            ++notCount;
            int notify = -2;
            MPI_Send(&notify, 1, MPI_INT, senderRank, toRequester, MPI_COMM_WORLD);

            if (notCount == size-1) {
                break;
            }
        }
        
        pthread_mutex_lock(&mutex);
        int checkSize = tasks.size();
        pthread_mutex_unlock(&mutex);
        int newTask;

        if (checkSize > recVar) {
            pthread_mutex_lock(&mutex);
            newTask = tasks.front();
            tasks.pop();
            pthread_mutex_unlock(&mutex);
        } else {
            if (checkSize < BORDER) {
                newTask = -2;
                ++notCount;
                end = true;
            } else {
                newTask = -1;
            }
        }

        MPI_Send(&newTask, 1, MPI_INT, senderRank, toRequester, MPI_COMM_WORLD);
        printf("Rank: %d -- send %d to %d\n", rank, newTask, senderRank);
    } while (true);
}

int main(int argc, char** argv) {
    int provided;

    MPI_Init_thread(&argc, &argv,  MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Error: MPI_THREAD_MULTIPLE not supported\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    double gloStart = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double* times;

    if (rank == 0) {
        times = new double[size];
    }

    int iterCounter = 0;
    pthread_attr_t attr;

    if (pthread_attr_init(&attr) != 0) {
        printf("Error: pthread_attr_init failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) != 0) {
        printf("Error: pthread_attr_setdetachstate failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    while (true) {
        for (int i = rank*100; i < (rank+1)*100; ++i) {
            int first = 50-i%100;
            if (first < 0) {
                first *= -1;
            }
            int second = rank-(iterCounter%size);
            if (second < 0) {
                second *= -1;
            }
            tasks.push(first*second*L);
        }

        othHave = true; iHave = true;
        int taskCounter = 0;

        if (pthread_create(&reqThread, &attr, request_tasks, NULL) != 0) {
            printf("Error: pthread_create failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (pthread_create(&handleThread, &attr, handle_requests, NULL) != 0) {
            printf("Error: pthread_create failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        printf("Begin of work\n");

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        while (othHave || iHave) {
            pthread_mutex_lock(&mutex);
            bool empty = tasks.empty();
            pthread_mutex_unlock(&mutex);

            if (empty) {
                iHave = false;
                //printf("Rank: %d -- empty\n", rank);
                continue;
            }

            pthread_mutex_lock(&mutex);
            int iter = tasks.front();
            tasks.pop();
            pthread_mutex_unlock(&mutex);

            taskCounter += iter;

            for (int i = 0; i < iter; ++i) {
                procRes += sin(i);
            }

            //printf("Rank: %d -- end of task\n", rank);
        }

        printf("Rank: %d -- end of work\n", rank);

        double end = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);

        if (pthread_join(reqThread, NULL) != 0) {
            printf("Error: pthread_join failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (pthread_join(handleThread, NULL) != 0) {
            printf("Error: pthread_join failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        printf("End of work\n");

        ++iterCounter;

        double diff = end-start;
        MPI_Gather(&diff, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            double max = *std::max_element(times, times+size);
            double disBalance = max - *std::min_element(times, times+size);

            printf("Iteration: %d, Disbalance: %f\n", iterCounter, disBalance);
            printf("Percentage of disbalance: %f\n", disBalance/max*100.0);
        }

        printf("Rank: %d, Iteration: %d, Time: %f, Result: %f, TaskCounter: %d\n",
            rank, iterCounter, end-start, procRes, taskCounter);

        if (iterCounter == 3) {
            break;
        }
    }

    if (rank == 0) {
        delete[] times;
    }

    double gloEnd = MPI_Wtime();

    if (rank == 0) {
        printf("Global time: %f\n", gloEnd-gloStart);
    }

    pthread_mutex_destroy(&mutex);
    MPI_Finalize();

    return EXIT_SUCCESS;
}