// caculate PI using MPI
// use the result of the Basel problem to find pi^2/6, and extract pi from it 

#include <stdio.h>
#include <mpi.h>
#include <math.h>
//#include <limits.h>

// take a small function to calculate N and stuff first, based on accuracy required
// let's say accuracy is 1E-14, find N such that 1/N^2 < 1E-14 

int main (int argc, char* argv[])
{
	double N=1E4;
	double accuracy;
	
    int rank, size, error;
	long long int i;
    double pi=0.0, result=0.0, sum=0.0, begin=0.0, end=0.0, delta=1E-14;
    
    error=MPI_Init (&argc, &argv);
    
    //Get process info
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
	
	if(rank==0)
	{
		scanf("%lf", &accuracy);
		N = 1.0/accuracy;
		if( N< 1E4 )
		{
			N=1E4;
		}
	}
    
    
    //Synchronize all processes and get the begin time
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
	
    //Each process will caculate a part of the sum
    for (i=rank; i<N; i+=size)
    {
		result+=1.0/( (i+1)*(i+1) );
    }
    
    //Sum up all results
    MPI_Reduce(&result, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //Synchronize all processes and get the end time
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    
    //Caculate and print pi
    if (rank==0)
    {
        pi=sqrt(sum*6);
        printf("np=%2d;    Time=%fs;    PI=%lf, N=%lf\n", size, end-begin, pi, N);
    }
    
    error=MPI_Finalize();
    
    return 0;
}