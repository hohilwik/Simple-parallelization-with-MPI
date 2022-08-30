#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include "matrix_multiply.h"

long **init_long_matrix(int rows, int cols)
{
  long *data = (long *) calloc(rows*cols, sizeof(long));
  long **array= (long **) calloc(rows, sizeof(long*));
  if(data == NULL || array == NULL) {
    fprintf(stderr, "Could not allocate sufficient memory!");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < rows; ++i) {
    array[i] = &(data[cols*i]);
  }
  return array;
}

void fprintf_matrix(FILE *stream, long** matrix, int rows, int cols)
{
  for(int i = 0; i < rows; ++i) 
  {
    for(int j = 0; j < cols; ++j) 
	{
        fprintf(stream, "%4ld", matrix[i][j]);
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}

int main(int argc, char *argv[])
{
  /* Get square matrix size */
  int N = DEFAULT_SIZE;

  /* a and b are the matrices to be multiplied, c is the result matrix. */
  long **a = NULL, **b = NULL, **c = NULL;

  // not using this anymore but the code breaks if I remove it?
  // so just letting it chill around
  int mul = 2;

  //variables used calculate elapsed time
  double t_start = 0, t_end = 0;
#ifdef BARRIERS
  double p_start = 0, p_end = 0;
#endif

  /* MPI related variables. */
  MPI_Status status;
  int my_rank, nproc;

  /* Auxiliary variables. */
  int averow, count, rows, extra, offset;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);

  /* Get MPI info. */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  

  /* Initialize matrices. */
  a = init_long_matrix(N, N);
  b = init_long_matrix(N, N);
  c = init_long_matrix(N, N);

  /* Initializing data. */
  if(my_rank == MASTER) {
	  
	printf("\nNumber of processes: %d\n", nproc);
    printf("Square matrix size: %d\n", N);
	
	srand(0);
    
    for(int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        a[i][j] = rand();
        b[i][j] = rand();
      }
    }

    printf("MPI_Wtime() precision: %lf\n", MPI_Wtick());
    t_start = MPI_Wtime();
  } 

  /* End of init step. Using synchronization barrier. */
#ifdef BARRIERS
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == MASTER) {
    p_start = MPI_Wtime();
  }
#endif

  /* Sending work to slaves. */
  if(my_rank == MASTER) {
    averow = N / nproc;
    extra = N % nproc;
    /* skipping the transfer between master and master cause that's silly. */
    offset = averow;
    for(int dest = 1; dest < nproc; ++dest) {
      rows = (dest <= extra) ? averow + 1 : averow;
      
      MPI_Send(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
      MPI_Send(&a[offset][0], rows*N, MPI_LONG, dest, FROM_MASTER,
               MPI_COMM_WORLD);
      MPI_Send(&(b[0][0]), N*N, MPI_LONG, dest, FROM_MASTER, MPI_COMM_WORLD);
      offset = offset + rows;
    }
    /* Afterwards, master should have rows=averows */
    rows = averow;
  } else {
      /* Receiving work from master. */
      MPI_Recv(&offset, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&rows, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      MPI_Recv(&a[0][0], rows*N, MPI_LONG, MASTER, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Get_count(&status, MPI_LONG, &count);
      
      MPI_Recv(&(b[0][0]), N*N, MPI_LONG, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_LONG, &count);
  }
  /* End of communication step. Using synchronization barrier. */
#ifdef BARRIERS
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == MASTER) {
    p_end = MPI_Wtime();
    printf("Communication step 1 elapsed time: %lf\n", p_end - p_start);
    p_start = MPI_Wtime();
  }
#endif


  /* Compute matrix multiplication */
  for(int k = 0; k < N; ++k) {
    for(int i = 0; i < rows; ++i) {
      for(int j = 0; j < N; ++j) {
        c[i][k] += a[i][j] * b[j][k];
      }
    }
  }

  /* End of computation step. Using synchronization barrier. */
#ifdef BARRIERS
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == MASTER) {
    p_end = MPI_Wtime();
    printf("Computation step elapsed time: %lf\n", p_end - p_start);
    p_start = MPI_Wtime();
  }
#endif

  /* Send results to master. */
  if(my_rank != MASTER) {
      MPI_Send(&offset, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
      MPI_Send(&c[0][0], rows*N, MPI_LONG, MASTER, FROM_WORKER, MPI_COMM_WORLD);
  }

  /* Receive results from worker tasks. */
  if (my_rank == MASTER) {
      for (int source = 1; source < nproc; ++source) {
        MPI_Recv(&offset, 1, MPI_INT, source, FROM_WORKER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rows, 1, MPI_INT, source, FROM_WORKER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&c[offset][0], rows*N, MPI_LONG, source, FROM_WORKER, 
          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  /* End of communication step. Using synchronization barrier. */
#ifdef BARRIERS
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == MASTER) {
    p_end = MPI_Wtime();
    printf("Communication step 2 elapsed time: %lf\n", p_end - p_start);
    p_start = MPI_Wtime();
  }
#endif

  if(my_rank == MASTER) {
    t_end = MPI_Wtime();
    printf("Total elapsed time: %lf\n", t_end - t_start);
  }


  MPI_Finalize();
  exit(EXIT_SUCCESS);
  //fancy return 0 ooof
}