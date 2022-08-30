#include<stdio.h>
#include<unistd.h>
#include<sys/syscall.h>
#include<stdlib.h>

int main()
{
	system("mpicc -o pi parallel_pi.c -lm -ldl");
	system("mpirun -np 1 ./pi");
	system("mpirun -np 2 ./pi");
	system("mpirun -np 3 ./pi");
	system("mpirun -np 4 ./pi");
	system("mpirun -np 5 ./pi");
	system("mpirun -np 6 ./pi");
	system("mpirun -np 7 ./pi");
	system("mpirun -np 8 ./pi");
	return 0;
	
}