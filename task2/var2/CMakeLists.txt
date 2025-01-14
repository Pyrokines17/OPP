cmake_minimum_required(VERSION 3.6)
project(var2 C)

set(CMAKE_C_STANDARD 17)

find_package(MPI REQUIRED)

add_executable(var2 main.c algebra.c)

target_link_libraries(var2 PRIVATE MPI::MPI_C)
target_link_libraries(var2 PRIVATE m)

include_directories(/usr/include/mpi/)