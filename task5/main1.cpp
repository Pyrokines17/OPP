#include <iostream>
#include <vector>
#include <mpi.h>

#define XLine 16
#define YLine 16

void initWorld(char*& world) {
    world[1] = 1;
    world[XLine+2] = 1;
    world[2*XLine] = 1;
    world[2*XLine+1] = 1;
    world[2*XLine+2] = 1;
}

void readPointNeighbors(int x, int y, std::vector<int>& neigh) {
    unsigned int k = 0;

    for (int i = x-1; i <= x+1; ++i) {
        for (int j = y-1; j <= y+1; ++j) {
            if (i == x && j == y) {
                continue;
            }

            if (i < 0) {
                neigh[2*k] = XLine+i;
            } else {
                neigh[2*k] = i%XLine;
            }

			neigh[2*k+1] = j;

            ++k;
        }
    }
}

unsigned int getLiveCountNeighbors(int x, int y, const char* world) {
    unsigned int count = 0;
    std::vector<int> neigh(8*2);
    int x1, y1;

    readPointNeighbors(x, y, neigh);

    for (int i = 0; i < 8; ++i) {
        x1 = neigh[i*2];
        y1 = neigh[i*2+1];

        if (world[y1*XLine+x1]) {
            ++count;
        }
    }

    return count;
}

void nextGenOfPart(char*& world, char* prevWorld, int size) {
    unsigned int ln;
    char p;

    for (int i = 2; i < YLine/size; ++i) {
        for (int j = 0; j < XLine; ++j) {
            p = prevWorld[i*XLine+j];
            ln = getLiveCountNeighbors(j, i, prevWorld);

            if (p) {
                if (ln < 2 || ln > 3) {
                    world[i*XLine+j] = 0;
                }
            } else {
                if (ln == 3) {
                    world[i*XLine+j] = 1;
                }
            }
        }
    }
}

void nextGenOfUp(char*& world, char* prevWorld) {
	unsigned int ln;
    char p;

    for (int i = 1; i < 2; ++i) {
        for (int j = 0; j < XLine; ++j) {
            p = prevWorld[i*XLine+j];
            ln = getLiveCountNeighbors(j, i, prevWorld);

            if (p) {
                if (ln < 2 || ln > 3) {
                    world[i*XLine+j] = 0;
                }
            } else {
                if (ln == 3) {
                    world[i*XLine+j] = 1;
                }
            }
        }
    }
}

void nextGenOfDown(char*& world, char* prevWorld, int size) {
	unsigned int ln;
    char p;

    for (int i = YLine/size; i < YLine/size+1; ++i) {
        for (int j = 0; j < XLine; ++j) {
            p = prevWorld[i*XLine+j];
            ln = getLiveCountNeighbors(j, i, prevWorld);

            if (p) {
                if (ln < 2 || ln > 3) {
                    world[i*XLine+j] = 0;
                }
            } else {
                if (ln == 3) {
                    world[i*XLine+j] = 1;
                }
            }
        }
    }
}

void cpPartOfWorld(char*& dest, const char* src, int size) {
    for (int i = XLine; i < XLine*(YLine/size+1); ++i) {
        dest[i] = src[i];
    }
}

void cpUpOfWorld(char*& dest, const char* src) {
    for (int i = 0; i < XLine; ++i) {
        dest[i] = src[i];
    }
}

void cpDownOfWorld(char*& dest, const char* src, int size) {
    for (int i = XLine*(YLine/size+1); i < XLine*(YLine/size+2); ++i) {
        dest[i] = src[i];
    }
}

void fullCP(char*& dest, const char* src, int size) {
    for (int i = 0; i < size; ++i) {
        dest[i] = src[i];
    }
}

char cmpPart(std::vector<char*>& w1, const char* w2, int size) {
    for (char* world : w1) {
        bool flag;

        for (int i = 0; i < size; ++i) {
            if (world[i] != w2[i]) {
                flag = false;
                break;
            } else {
                flag = true;
                continue;
            }
        }

        if (flag) {
            return 1;
        } else {
            continue;
        }
    }

    return 0;
}

void printWorld(const char* world) {
    for (int i = 0; i < XLine*2; ++i) {
        std::cout << "-";
    } std:: cout << std::endl;

    for (int i = 0; i < YLine; ++i) {
        for (int j = 0; j < XLine; ++j) {
            if (world[i*XLine+j]) {
                std::cout << "@";
            } else {
                std::cout << "*";
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < XLine*2; ++i) {
        std::cout << "-";
    } std:: cout << std::endl;
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int size, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char flag = 1, locOpt = 0, opt = 0;
	char* world,* partOfWorld,* prevPartOfWorld;

	partOfWorld = new char[XLine*(YLine/size+2)]{0};
	prevPartOfWorld = new char[XLine*(YLine/size+2)]{0};

	if (rank == 0) {
		world = new char[XLine*YLine]{0};
		initWorld(world);
		printWorld(world);
	}

	MPI_Scatter(world, XLine*YLine/size, MPI_CHAR, 
		&partOfWorld[XLine], XLine*YLine/size, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Request SReq, SReq1;
	MPI_Request RReq, RReq1;

	std::vector<char*> parts;

	while (flag) {
		MPI_Isend(&partOfWorld[XLine], XLine, MPI_CHAR, (rank+size-1)%size, 0, MPI_COMM_WORLD, &SReq);
		MPI_Isend(&partOfWorld[XLine*YLine/size], XLine, MPI_CHAR, (rank+1)%size, 1, MPI_COMM_WORLD, &SReq1);

		MPI_Irecv(&partOfWorld[0], XLine, MPI_CHAR, (rank+size-1)%size, 1, MPI_COMM_WORLD, &RReq);
		MPI_Irecv(&partOfWorld[XLine*(YLine/size+1)], XLine, MPI_CHAR, (rank+1)%size, 0, MPI_COMM_WORLD, &RReq1);

		cpPartOfWorld(prevPartOfWorld, partOfWorld, size);
		nextGenOfPart(partOfWorld, prevPartOfWorld, size);

		MPI_Wait(&SReq, MPI_STATUS_IGNORE);
		MPI_Wait(&RReq, MPI_STATUS_IGNORE);

		cpUpOfWorld(prevPartOfWorld, partOfWorld);
		nextGenOfUp(partOfWorld, prevPartOfWorld);

		MPI_Wait(&SReq1, MPI_STATUS_IGNORE);
		MPI_Wait(&RReq1, MPI_STATUS_IGNORE);

		cpDownOfWorld(prevPartOfWorld, partOfWorld, size);
		nextGenOfDown(partOfWorld, prevPartOfWorld, size);

		MPI_Gather(&partOfWorld[XLine], XLine*YLine/size, MPI_CHAR, world, XLine*YLine/size, MPI_CHAR, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printWorld(world);
		}

		char* part1 = new char[XLine*YLine/size];
		fullCP(part1, prevPartOfWorld+XLine, XLine*YLine/size);

		parts.push_back(part1);
		locOpt = cmpPart(parts, partOfWorld+XLine, XLine*YLine/size);
		MPI_Allreduce(&locOpt, &opt, 1, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD);

		if (opt == size) {
			flag = 0;
		} else {
			opt = 0;
		}
	}

	for (char* obj : parts) {
        delete[] obj;
    }

	delete[] partOfWorld;
	delete[] prevPartOfWorld;

	if (rank == 0) {
		delete[] world;
	}

	MPI_Finalize();
	return 0;
}