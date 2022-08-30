#include<stdio.h>
#include<unistd.h>
#include<sys/syscall.h>
#include<stdlib.h>

int main()
{
	//system("mpicc -o pi parallel_pi.c -lm -ldl");
	//system("mpirun -np 1 ./matrix_multiply");
	//system("mpirun -np 2 ./matrix_multiply");
	//system("mpirun -np 3 ./matrix_multiply");
	//system("mpirun -np 4 ./matrix_multiply");
	//system("mpirun -np 5 ./matrix_multiply");
	system("mpirun -np 6 ./matrix_multiply");
	system("mpirun -np 7 ./matrix_multiply");
	system("mpirun -np 8 ./matrix_multiply");
	return 0;
	
}