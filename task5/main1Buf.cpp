#include <iostream>
#include <vector>
#include <mpi.h>

#define XLine 250
#define YLine 250

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

    for (int i = 2; i < size; ++i) {
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

    for (int i = size; i < size+1; ++i) {
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
    for (int i = XLine; i < XLine*(size+1); ++i) {
        dest[i] = src[i];
    }
}

void cpUpOfWorld(char*& dest, const char* src) {
    for (int i = 0; i < XLine; ++i) {
        dest[i] = src[i];
    }
}

void cpDownOfWorld(char*& dest, const char* src, int size) {
    for (int i = XLine*(size+1); i < XLine*(size+2); ++i) {
        dest[i] = src[i];
    }
}

void fullCP(char*& dest, const char* src, int size) {
    for (int i = 0; i < size; ++i) {
        dest[i] = src[i];
    }
}

void cmpPart(std::vector<char*>& w1, const char* w2, int size, char*& locFlags) {
    locFlags[0] = 0;

    for (int i = 1; i <= w1.size(); ++i) {
        bool flag;

        for (int j = 0; j < size; ++j) {
            if (w1[i-1][j] != w2[j]) {
                flag = false;
                break;
            } else {
                flag = true;
                continue;
            }
        }

        if (flag) {
            locFlags[i] = 1;
        } else {
            locFlags[i] = 0;
        }
    }
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

    if (YLine < size) {
        if (rank == 0) {
            std::cout << "Need more string's (minStr == sizeOfThr)" << std::endl;
        }

        MPI_Finalize();
        return 0;
    }

	char flag = 1, locOpt = 0, opt = 0;
	char* world = nullptr,* partOfWorld = nullptr,* prevPartOfWorld = nullptr;

    int sizeOfPart = (rank+1) * YLine/size - rank * YLine/size;
    //std::cout << sizeOfPart << ":" << rank << std::endl;

	partOfWorld = new char[XLine*(sizeOfPart+2)]{0};
	prevPartOfWorld = new char[XLine*(sizeOfPart+2)]{0};

	if (rank == 0) {
		world = new char[XLine*YLine]{0};
		initWorld(world);
		
        //printWorld(world);
	}

    int* begins = nullptr;
    int* counts = nullptr;

    if (rank == 0) {
        begins = new int[size]{0};
        counts = new int[size]{0};

        for (int i = 0; i < size; ++i) {
            begins[i] = i * YLine / size * XLine;
            counts[i] =(i+1) * YLine/size * XLine - begins[i];
        } 
    }

	//MPI_Scatter(world, XLine*YLine/size, MPI_CHAR, &partOfWorld[XLine], XLine*YLine/size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatterv(world, counts, begins, MPI_CHAR, &partOfWorld[XLine], sizeOfPart*XLine, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Request SReq, SReq1;
	MPI_Request RReq, RReq1;
    MPI_Request red;

	std::vector<char*> parts;

    int counter = 1;
    double begin = MPI_Wtime();

    char* locFlags = (char*)calloc(1, sizeof(char));
    char* gloFlags = (char*)calloc(1, sizeof(char));

    int dims[1] = {size}, periods[1] = {1};
    int left, right;
    MPI_Comm ring;

    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &ring);
    MPI_Cart_shift(ring, 0, 1, &left, &right);

	while (flag) {
        MPI_Iallreduce(locFlags, gloFlags, counter, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD, &red);

		MPI_Isend(&partOfWorld[XLine], XLine, MPI_CHAR, left, 0, ring, &SReq);
		MPI_Isend(&partOfWorld[XLine*sizeOfPart], XLine, MPI_CHAR, right, 1, ring, &SReq1);

		MPI_Irecv(&partOfWorld[0], XLine, MPI_CHAR, left, 1, ring, &RReq);
		MPI_Irecv(&partOfWorld[XLine*(sizeOfPart+1)], XLine, MPI_CHAR, right, 0, ring, &RReq1);

		cpPartOfWorld(prevPartOfWorld, partOfWorld, sizeOfPart);
		nextGenOfPart(partOfWorld, prevPartOfWorld, sizeOfPart);

		MPI_Wait(&SReq, MPI_STATUS_IGNORE);
		MPI_Wait(&RReq, MPI_STATUS_IGNORE);

		cpUpOfWorld(prevPartOfWorld, partOfWorld);
		nextGenOfUp(partOfWorld, prevPartOfWorld);

		MPI_Wait(&SReq1, MPI_STATUS_IGNORE);
		MPI_Wait(&RReq1, MPI_STATUS_IGNORE);

		cpDownOfWorld(prevPartOfWorld, partOfWorld, sizeOfPart);
		nextGenOfDown(partOfWorld, prevPartOfWorld, sizeOfPart);

        /*
		//MPI_Gather(&partOfWorld[XLine], XLine*sizeOfPart, MPI_CHAR, world, XLine*sizeOfPart, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&partOfWorld[XLine], sizeOfPart*XLine, MPI_CHAR, world, counts, begins, MPI_CHAR, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printWorld(world);
		}
        */

		char* part1 = new char[XLine*sizeOfPart];
		fullCP(part1, prevPartOfWorld+XLine, XLine*sizeOfPart);
		parts.push_back(part1);
        
        MPI_Wait(&red, MPI_STATUS_IGNORE);
        
        for (int i = 0; i < counter; ++i) {
            if (gloFlags[i] == size) {
                flag = 0;

                if (rank == 0) {
                    std::cout << 
                        "Optimal configuration detected" 
                            << std::endl;
                }

                break;
            } else continue;
        }

        ++counter;
        free(locFlags); free(gloFlags);

        gloFlags = (char*)calloc(counter, sizeof(char));
        locFlags = (char*)malloc(sizeof(char) * counter);

		cmpPart(parts, partOfWorld+XLine, XLine*sizeOfPart, locFlags);
	}

    free(locFlags); free(gloFlags);
    double end = MPI_Wtime();

	for (char* obj : parts) {
        delete[] obj;
    }

	delete[] partOfWorld;
	delete[] prevPartOfWorld;

	if (rank == 0) {
		delete[] world;
        delete[] begins;
        delete[] counts;

        std::cout << end-begin << std::endl;
	}

	MPI_Finalize();
	return 0;
}