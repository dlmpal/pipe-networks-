cmake_minimum_required(VERSION 3.21)
project(PipeNetwork C)

set(CMAKE_C_STANDARD 99)

add_executable(PipeNetwork main.c matrix_assembly.c matrix_assembly.h linalg.c linalg.h)

