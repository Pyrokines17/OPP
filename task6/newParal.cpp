#include <algorithm>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <atomic>
#include <math.h>
#include <mpi.h>
#include <queue>

#define L 200
#define MIN 25
#define BORDER 5
#define MAX_ITER 5

int rank, size;
double procRes = 0;

std::queue<int> tasks;
pthread_mutex_t mutex;

std::atomic<bool> startReq, startHandle,
    needReq, needHandle, isEmpty;

pthread_t reqThread, handleThread;

const int toHandler = 0, toRequester = 1;

void* request_tasks(void* arg) {
    while (needReq) {
        if (!startReq) {
            continue;
        }

        if (size == 1) {
            startReq = false;
            continue;
        }

        bool* havent = new bool[size]{false};
        int countHavent = 0;

        while (true) {
            pthread_mutex_lock(&mutex);
            int mySize = tasks.size();
            pthread_mutex_unlock(&mutex);

            if (mySize > MIN) {
                continue;
            }

            for (int i = 1; i < size; ++i) {
                int nextRank = (rank+i)%size;
                int newTask;

                if (havent[nextRank]) {
                    continue;
                }

                MPI_Send(&mySize, 1, MPI_INT, nextRank, toHandler, MPI_COMM_WORLD);
                MPI_Recv(&newTask, 1, MPI_INT, nextRank, toRequester, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (newTask == -2) {
                    havent[nextRank] = true;
                    ++countHavent;
                    continue;
                } else if (newTask == -1) {
                    continue;
                } else {
                    pthread_mutex_lock(&mutex);
                    tasks.push(newTask);
                    pthread_mutex_unlock(&mutex);
                    isEmpty = false;
                    continue;
                }
            }

            if (countHavent == size-1) {
                startReq = false;
                break;
            }
        }
    }

    return EXIT_SUCCESS;
}

void* handle_requests(void* arg) {
    while (needHandle) {
        if (!startHandle) {
            continue;
        }

        if (size == 1) {
            startHandle = false;
            continue;
        }

        int recvSize;
        bool end = false;

        MPI_Status status;

        int notify = -2;
        int notCount = 0;

        while (true) {
            MPI_Recv(&recvSize, 1, MPI_INT, MPI_ANY_SOURCE, toHandler, MPI_COMM_WORLD, &status);
            int senderRank = status.MPI_SOURCE;

            if (end) {
                MPI_Send(&notify, 1, MPI_INT, senderRank, toRequester, MPI_COMM_WORLD);
                ++notCount;

                if (notCount == size-1) {
                    startHandle = false;
                    break;
                }

                continue;
            }

            pthread_mutex_lock(&mutex);
            int mySize = tasks.size();
            pthread_mutex_unlock(&mutex);
            int task;

            if (mySize > recvSize) {
                pthread_mutex_lock(&mutex);
                task = tasks.front();
                tasks.pop();
                pthread_mutex_unlock(&mutex);
            } else {
                if (mySize < BORDER) {
                    end = true;
                } 
                
                task = -1;
            }

            MPI_Send(&task, 1, MPI_INT, senderRank, toRequester, MPI_COMM_WORLD);
        }
    }

    return EXIT_SUCCESS;
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
    pthread_mutexattr_t mutexAttr;

    if (pthread_attr_init(&attr) != 0) {
        printf("Error: pthread_attr_init\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (pthread_mutexattr_init(&mutexAttr) != 0) {
        printf("Error: pthread_mutexattr_init\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    needReq = true;
    startReq = false;

    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) != 0) {
        printf("Error: pthread_attr_setdetachstate\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (pthread_create(&reqThread, &attr, request_tasks, NULL) != 0) {
        printf("Error: pthread_create\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    needHandle = true;
    startHandle = false;

    if (pthread_create(&handleThread, &attr, handle_requests, NULL) != 0) {
        printf("Error: pthread_create\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (pthread_mutex_init(&mutex, &mutexAttr) != 0) {
        printf("Error: pthread_mutex_init\n");
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

        int taskCounter = 0;
        startReq = true;
        startHandle = true;
        isEmpty = false;

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        while (startReq || startHandle || !isEmpty) {
            pthread_mutex_lock(&mutex);
            int mySize = tasks.size();
            pthread_mutex_unlock(&mutex);

            if (mySize == 0) {
                isEmpty = true;
                continue;
            }

            pthread_mutex_lock(&mutex);
            int task = tasks.front();
            tasks.pop();
            pthread_mutex_unlock(&mutex);

            taskCounter += task;
            
            for (int i = 0; i < task; ++i) {
                procRes += sin(i);
            }
        }

        double end = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);

        ++iterCounter;

        double diff = end-start;
        MPI_Gather(&diff, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            double max = *std::max_element(times, times+size);
            double disBalance = max - *std::min_element(times, times+size);

            printf("Iteration: %d, max: %f, disbalance: %f\n", iterCounter, max, disBalance);
            printf("Percentage of disbalance: %f\n", disBalance/max*100.0);
        }

        printf("Rank: %d, Iteration: %d, Time: %f, Result: %f, TaskCounter: %d\n", 
            rank, iterCounter, end-start, procRes, taskCounter);

        if (iterCounter == MAX_ITER) {
            needReq = false;
            needHandle = false;
            break;
        }
    }

    if (pthread_join(reqThread, NULL) != 0) {
        printf("Error: pthread_join failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (pthread_join(handleThread, NULL) != 0) {
        printf("Error: pthread_join failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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