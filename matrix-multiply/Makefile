CC = mpicc
CFLAGS = -std=gnu99 -O3 -Wall
# OMPI_MPICC = clang

all: matrix_multiply

matrix_multiply: parallel_matrix_multiply.c matrix_multiply.h
ifdef OMPI_MPICC
	env OMPI_MPICC=$(OMPI_MPICC) $(CC) $(CFLAGS) -o matrix_multiply parallel_matrix_multiply.c
else
	$(CC) $(CFLAGS) -o matrix_multiply parallel_matrix_multiply.c
endif

clean:
	rm -f matrix_multiply