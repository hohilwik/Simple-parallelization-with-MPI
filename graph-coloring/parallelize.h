#ifndef PARALLELIZE
#define PARALLELIZE

#define TAG 42
#define POS(x, y, side) ((x) + (y) * (side))
#include<mpi.h>
#include<stdio.h>
#include "hashset.h"

int get_nodes_per_process(int nodes, int max_processes);
int get_pid_of_node(int node, int nodes, int max_processes);
int rand_between(int low, int high);
int color_node(int node, int pid, int nodes_per_block, int * madjacency, int* colors, CustomSet * psetcolors, int nodes);
int is_local(int node, int i, int * madjacency, int nodes, int p);


int plassman_algorithm_p_processors(int argc, char * argv[]);
void plassman_p_processors(int pid, int p, MPI_Comm * comm_workers);
void plassman_v_processors(int pid);
int find_my_color(CustomSet * psetcolors, int nodes);
void print_adjacency_matrix(int pid, int * madjacency, int my_nodes, int nodes) ;
void error(char * msg);
void read_graph_to_adjacency_matrix(FILE * f, int ** pmatrix, int * pnodes);

#endif