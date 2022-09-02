# Simple-parallelization-with-MPI
Written as part of a Distributed Systems assignment. Distributed algorithms for calculating pi, multiplying matrices, and graph coloring. Uses Jones-Plassman for the graph coloring problem. Deleted previous Cole-Vishkin code.

## Q2

Calculating pi from the result of the Basel problem:

- Number of terms was decided by desired accuracy
- For N=10^5 terms, the context switching time between different processes far outweighs the compute benefit of multiple processes. Time taken goes from 0.0012s to 0.061s for two processes, and then goes down to 0.015s for 8 processes.
- For N=10^7 terms, the compute power benefit of multiple processors is clearly seen, time taken going from 0.1214s for 1 process to 0.0154s for 8 processes

### How to run

- ``cd ./calc-pi``
- ``gcc run.c -o run``
- ``./run`` 


## Q4

Distributed Matrix multiplication

- For a small matrix of N=40, multiple processes slow down the computation, going from 0.000268s to 0.000427s from 1 process to 2 processes.
- For larger matrices, N=200, time taken goes from 0.035s for 1 process to 0.007s for 8 processes
- For much larger matrices, N=4096, time taken goes from 1168s for 1 process, to 193s for 8 processes.

### How to run

- ``cd ./matrix-multiply``
- ``make``
- ``gcc run.c -o run``
- ``./run``
- change DEFAULT_SIZE in `matrix_multiply.h` to change matrix size

## Q6

Distributed graph coloring

- For a small planar graph with 50000 nodes, it uses 5 colors, close to the optimal 3, and the ideal 3+1 for the algorithm
- For graphs with less than 100,000 nodes, speed up for 4 cores vs 1 core is not seen
- For graphs with 10 million nodes, a nearly linear speed up is noticed, with time taken going from 64954s for 1 thread to 7824s for 8 threads

### How to run

- ``cd ./graph-coloring``
- ``python3 gen_graphs.py <size> <probability of nodes being connected in range(0,1)> > graph.txt
- ``gcc run.c -o run``
- ``./run``
- use graph_gen.py to generate graphs of different sizes



