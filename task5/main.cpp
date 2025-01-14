#include <iostream>
#include <vector>

#define XLine 10
#define YLine 10

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

            if (j < 0) {
                neigh[2*k+1] = YLine+j;
            } else {
                neigh[2*k+1] = j%YLine;
            }

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

void nextGen(char*& world, char* prevWorld) {
    unsigned int ln;
    char p;

    for (int i = 0; i < YLine; ++i) {
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

void cpWorld(char*& dest, const char* src) {
    for (int i = 0; i < XLine*YLine; ++i) {
        dest[i] = src[i];
    }
}

bool cmpWorld(std::vector<char*>& w1, const char* w2) {
    for (char* world : w1) {
        bool flag;

        for (int i = 0; i < XLine*YLine; ++i) {
            if (world[i] != w2[i]) {
                flag = false;
                break;
            } else {
                flag = true;
                continue;
            }
        }

        if (flag) {
            return true;
        } else {
            continue;
        }
    }

    return false;
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

int main() {
    char* world = new char[XLine*YLine]{0};
    char* prevWorld = new char[XLine*YLine]{0};

    initWorld(world);

    bool isOpt = false;

    std::vector<char*> worlds;
    
    //printWorld(world);

    while (!isOpt) {
        cpWorld(prevWorld, world);
        nextGen(world, prevWorld);

        //printWorld(world);

        char* world1 = new char[XLine*YLine]{0};
        cpWorld(world1, prevWorld);

        worlds.push_back(world1);
        isOpt = cmpWorld(worlds, world);

        if (isOpt) {
            std::cout << 
                "Optimal configuration detected" 
                    << std::endl;
        }
    }

    for (char* obj : worlds) {
        delete[] obj;
    }

    delete[] world;
    delete[] prevWorld;

    return 0;
}
