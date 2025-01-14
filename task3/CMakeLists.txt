cmake_minimum_required(VERSION 3.6)
project(var1 C)

set(CMAKE_C_STANDARD 17)

add_executable(var1 main.c algebra.c)

target_link_libraries(var1 PRIVATE m)