#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include "parallelize.h"
#include "hashset.h"

#define MIN(x,y) (x<y?x:y)

void read_graph_to_adjacency_matrix(FILE * f, int ** pmatrix, int * pnodes) {
	int edges, i;
	fscanf(f, " %d", pnodes);
	fscanf(f, " %d", &edges);
	*pmatrix = (int *) calloc(*pnodes * *pnodes, sizeof(int));
	if (*pmatrix == NULL) {
		error("Out of memory in read_graph_to_adjacency_matrix()\n");
		MPI_Finalize();
		exit(-1);
	}
	
	/* read edges one by one */
	for (i = 0; i < edges; i++) {
		int node1, node2;
		int n = fscanf(f, " %d %d", &node1, &node2);
		if (n != 2 || node1 < 0 || node1 >= *pnodes || node2 < 0 || node2 >= *pnodes) {
			fprintf(stderr, "Wrong edge %d in the graph definition file\n", i);
			exit(1);
		}
		(*pmatrix)[POS(node1, node2, *pnodes)] = 1;
		(*pmatrix)[POS(node2, node1, *pnodes)] = 1;
	}
}

int get_nodes_per_process(int nodes, int max_processes) {
	return floor( nodes * 1.0 / max_processes );
}

int get_pid_of_node(int node, int nodes, int max_processes) {
	return MIN(floor( node * 1.0 / get_nodes_per_process(nodes, max_processes ) ), max_processes-1) + 1;
}



int find_my_color(CustomSet * psetcolors, int nodes) {
	int color;
	for (color = 1; customSetContains(psetcolors, (void *) color) && color < nodes; color++);
	return color;
}

/*
 * Generate a random value between a and b inclusive. It suposes the srand function was already called.
 */
int rand_between(int low, int high) {

	return rand() % (high - low + 1) + low;
}


int color_node(int node, int pid, int nodes_per_block, int * madjacency, int* colors, CustomSet * psetcolors, int nodes) {
	int othernode;
	customSetClear(psetcolors);
	for (othernode = 0; othernode < nodes; othernode++) {
		if (othernode != (pid-1)*nodes_per_block + node) {
			if (madjacency[POS(othernode, node, nodes)] == 1) {
				if (colors[othernode] != 0)
				/* the neighbor already has a color */
					customSetAdd(psetcolors, (void *) colors[othernode]);
			}
		}
	}

	return find_my_color(psetcolors, nodes);
}


int is_local(int node, int i, int * madjacency, int nodes, int p) {
	int othernode;
	int pid = get_pid_of_node(node, nodes, p);
	for (othernode = 0; othernode < nodes; othernode++) {
		if (othernode != node && madjacency[POS(othernode, i, nodes)] == 1) {
			if (get_pid_of_node(othernode, nodes, p) != pid)
				return 0;
		}
	}
	return 1;
}


void print_adjacency_matrix(int pid, int * madjacency, int my_nodes, int nodes) {
	int i, j;
	printf("[rank %d] adjacency %d:\n", pid, my_nodes*nodes);
	for ( i = 0; i < my_nodes; i++ ) {
		for ( j = 0; j < nodes; j++ )	
			printf( "%d ", madjacency[i*nodes + j]);
		printf( "[rank %d]\n",pid);
	}
	
	fflush(stdout);
}

int plassman_algorithm_p_processors(int argc, char * argv[]) {
	int * matrix, * colors, * conflicts;
	int numnodes, numcolors;
	int i, j, err, blocks;
	int  numtasks, rank, conflict_pos = 0, read_items;
	static int excl_ranks[1] = {0};
	clock_t ticks1, ticks2;
	MPI_Request * reqs;
	MPI_Status * status;

	MPI_Group orig_group, workers_group;
	MPI_Comm comm_workers;

	CustomSet * psetcolors = (CustomSet *) customSetNewHashSet();
    ticks1 = clock();
	
	err = MPI_Init(&argc, &argv); /* Initialize MPI */
	if (err != MPI_SUCCESS) {
		printf("MPI_init failed!\n"); return 1;
	}

	err = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	/* Get nr of tasks */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_size failed!\n"); return 1;
	}

	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get id of this process */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_rank failed!\n"); return 1;
	}
	

	// create the workers group = all \ {master}
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	MPI_Group_excl(orig_group, 1, excl_ranks, &workers_group);

	MPI_Comm_create(MPI_COMM_WORLD, workers_group, &comm_workers);

	if (rank == 0) {
		read_graph_to_adjacency_matrix(stdin, &matrix, &numnodes);

		colors = (int*)malloc(numnodes*sizeof(int));
		conflicts = (int*)malloc(numnodes*sizeof(int));

		reqs = (MPI_Request*)malloc((numtasks-1)*sizeof(MPI_Request));
		status = (MPI_Status*)malloc((numtasks-1)*sizeof(MPI_Status));
        
		blocks = floor(numnodes *1.0/ (numtasks-1));
		for (i = 1; i < numtasks-1; i++){
			MPI_Send( &numnodes, 1, MPI_INT, i, TAG, MPI_COMM_WORLD );
			MPI_Send( &matrix[(i-1) * numnodes * blocks], numnodes*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD );
		}

		MPI_Send( &numnodes, 1, MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );
		MPI_Send( &matrix[(numtasks-2) * numnodes * blocks], numnodes*(numnodes-(numtasks-2)*blocks), MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );

		for (i = 1; i < numtasks-1; i++)
			MPI_Irecv( &colors[(i-1)*blocks], blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[i-1] );
		MPI_Irecv( &colors[(numtasks-2)*blocks], numnodes - (numtasks-2)*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[numtasks-2] );
		MPI_Waitall(numtasks-1, reqs, status);
		
	
		numcolors = 0;
		// FIND MAX COLOR
		for (i = 1; i < numnodes; i++)
			if (colors[i] > numcolors)   
				numcolors = colors[i];

		printf("%d\n%d\n", numnodes, numcolors);
		for (i = 0; i < numnodes; i++) {
			printf("%d\n", colors[i]);
		}
        ticks2	= clock();
        printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);

	} else {	
		plassman_p_processors(rank, numtasks-1, &comm_workers);
	}
	

	err = MPI_Finalize();	         /* Terminate MPI */
	if (err != MPI_SUCCESS) {
		printf("Error in MPI_Finalize!\n");
		return -1;
	}
	
    return 0;
}

void plassman_p_processors(int pid, int p, MPI_Comm * comm_workers) {
	int nodes;
	int * madjacency;
	int * all_weights, i, j;
	int my_nodes, nodes_per_block, max_nodes_per_block, n_locals;
	int * send_to, * n_wait, * n_send, * locals;
	int * receive_from, ended = 0;
	int * colors, my_colors[2], tmp = 0;

	MPI_Request req;
	MPI_Status status;
	CustomSet * psetcolors = (CustomSet *) customSetNewHashSet();
	int good = 1;	

	MPI_Recv( &nodes, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	nodes_per_block = get_nodes_per_process(nodes, p);
	max_nodes_per_block = nodes - (p - 1) * nodes_per_block;
	

	if (pid < p)
		my_nodes = nodes_per_block;
	else
		my_nodes = max_nodes_per_block;
		
	n_locals = 0;
	//printf("[rank %d] nodes: %d\n", pid, nodes);

	n_wait = (int*) malloc( my_nodes * sizeof(int) );
	n_send = (int*) malloc( my_nodes * sizeof(int) );
	send_to = (int*) malloc( my_nodes*(nodes-1) * sizeof(int) );
	receive_from = (int*) malloc( my_nodes*(nodes-1) * sizeof(int) );
	
	colors = (int*) malloc( (nodes) * sizeof(int) );
	all_weights = (int*) malloc( (nodes) * sizeof(int) );
	madjacency = (int*) malloc( nodes * my_nodes * sizeof(int) );
	locals = (int*) malloc( my_nodes * sizeof(int) );
	MPI_Recv( madjacency, nodes*my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	//random allocation
	srand( time( NULL ) / pid );
	for (i = 0; i < my_nodes; i++){
		all_weights[(pid-1)*nodes_per_block +i] = rand();
		n_wait[i] = 0;
		n_send[i] = 0;
	}
		

	//exchange all weights
	for ( i = 0; i < p-1; i++)
		MPI_Bcast( &all_weights[i*nodes_per_block], nodes_per_block, MPI_INT, i, *comm_workers );
	MPI_Bcast( &all_weights[(p-1)*nodes_per_block], max_nodes_per_block, MPI_INT, p-1, *comm_workers );
	
	//printf("[rank %d] weights exchanged\n", pid);

	//fflush(stdout);

	for ( i = 0; i < my_nodes; i++ ) {
		for (j = 0; j < nodes; j++) {

			if ((pid-1)*nodes_per_block + i == j || madjacency[i*nodes + j] <= 0 || pid == get_pid_of_node(j, nodes, p))
				continue;

			if ( all_weights[(pid-1)*nodes_per_block + i] < all_weights[j] ) {
				receive_from[i*(nodes-1) + n_wait[i]] = j;
				n_wait[i]++;
			} 
			else {
				send_to[i*(nodes-1) + n_send[i]] = j;
				n_send[i]++;
			}
		}

		locals[i] = is_local( (pid-1)*nodes_per_block+i, i, madjacency, nodes, p);
		if (locals[i]) {
			n_locals++;
			continue;
		}

		//printf("[rank %d] vertice[%d]: wait %d | send %d\n", pid, (pid-1)*nodes_per_block+i, n_wait[i], n_send[i]);
		if(n_wait[i] == 0 ){
			colors[(pid-1)*nodes_per_block+i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
			ended++;
			
			my_colors[0] = (pid-1)*nodes_per_block+i;
			my_colors[1] = colors[my_colors[0]];
			for (j = 0; j < n_send[i]; j++) {
				MPI_Isend(my_colors, 2, MPI_INT, get_pid_of_node(send_to[i*(nodes-1) + j], nodes, p), TAG, MPI_COMM_WORLD, &req);
			}
		}
	}
	
	//print_adjacency_matrix(pid, madjacency, my_nodes, nodes);

	while(ended != my_nodes-n_locals) {

		MPI_Recv(my_colors, 2, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);

		colors[my_colors[0]] = my_colors[1];

		for ( i = 0; i < my_nodes; i++ ) {
			if (n_wait[i] == 0 || locals[i])
				continue;
			//printf("[rank %d] %d not done...\n", pid, i);
			for (j = 0; j < n_wait[i]; j++)
				if (receive_from[i*(nodes-1)+j] == my_colors[0] ) {
					receive_from[i*(nodes-1)+j] = receive_from[i*(nodes-1)+n_wait[i]-1];
					n_wait[i]--;
					break;
				}
				
			if (n_wait[i] > 0)
				continue;
			//printf("[rank %d] %d done...\n", pid, i);
			colors[(pid-1)*nodes_per_block+i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
			ended++;
			
			my_colors[0] = (pid-1)*nodes_per_block+i;
			my_colors[1] = colors[my_colors[0]];
			for (j = 0; j < n_send[i]; j++) {
				MPI_Isend(my_colors, 2, MPI_INT, get_pid_of_node(send_to[i*(nodes-1) + j], nodes, p), TAG, MPI_COMM_WORLD, &req);
			}
		}

	}

	for(i = 0; i < my_nodes; i++)
		if (locals[i])
			colors[(pid-1)*nodes_per_block + i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
		//printf("[rank %d] color[%d] = %d\n", pid, i, colors[(pid-1)*nodes_per_block + i]);
	//printf("[rank %d] Send colors to master\n", pid);
	MPI_Send( &colors[(pid-1)*nodes_per_block], my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD);
}

int main(int argc, char * argv[]) 
{
    plassman_algorithm_p_processors(argc,argv);
    return 0;
}