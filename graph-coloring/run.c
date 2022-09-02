#include<stdio.h>
#include<unistd.h>
#include<sys/syscall.h>
#include<stdlib.h>

int main()
{
	system("mpicc -o graphcolor hashset.c hashset.h parallelize.h parallelize.c -lm");
	system("mpirun -np 1 ./graphcolor");
	system("mpirun -np 2 ./graphcolor");
	system("mpirun -np 3 ./graphcolor");
	system("mpirun -np 4 ./graphcolor");
	system("mpirun -np 5 ./graphcolor");
	system("mpirun -np 6 ./graphcolor");
	system("mpirun -np 7 ./graphcolor");
	system("mpirun -np 8 ./graphcolor");
	return 0;
	
}