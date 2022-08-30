#!/bin/sh
# run mpi_pi.c with np from 1 to 8

mpicc -o pi parallel_pi.c -lm -ldl


mpirun -np 8 ./pi
