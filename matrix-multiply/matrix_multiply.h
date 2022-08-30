#ifndef _MATRIX_MUL_H
 #define _MATRIX_MUL_H
  #define BARRIERS
  //#define DEBUG
  //#define VERBOSE
  #define MASTER 0
  #define FROM_MASTER 1
  #define FROM_WORKER 2
  #define DEFAULT_SIZE 2048
#endif

long **init_long_matrix(int rows, int cols);
void fprintf_matrix(FILE *stream, long** matrix, int rows, int cols);