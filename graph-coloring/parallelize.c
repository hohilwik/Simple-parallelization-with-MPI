#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

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

int main(int argc, char * argv[]) 
{
    plassman_algorithm_p_processors(argc,argv)
    return 0;
}