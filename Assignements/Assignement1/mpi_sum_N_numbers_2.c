#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define USE MPI

int main ( int argc , char *argv[ ] )
{

  unsigned long long int sum=0, local_sum ; 
   
  double start_time, end_time;   
  int myid , numprocs , proc ;
  MPI_Status status;
  MPI_Request request;
  // master process
  int master = 0;
  int tag = 123;
  unsigned long long int input;

  FILE* in_file = fopen("input_parallel_sum.txt", "r"); 
         
  if (! in_file )
    {  
      printf("Error reading file\n"); 
      exit(-1); 
    } 

  fscanf(in_file, "%llu", &input );
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);


  unsigned int* array=(unsigned  int*) calloc(input, sizeof(unsigned int));
  unsigned int i;
  unsigned int resto = input%numprocs;
  for (i = 0; i < input; i++)
    array[i] = i + 1;
  unsigned int N = input/numprocs;
// take time of processors after initial I/O operation
  start_time = MPI_Wtime();
  local_sum=0;
  for (i=0; i<N ; i++) {
    local_sum+=array[myid*N+i];
  }
  if (myid==0)
    {
      printf("\n # partial sum on master processor: %llu\n", local_sum);
    }
  else
    {
      printf ("\n # partial sum on processor %d: %llu\n", myid, local_sum);
	}
  MPI_Reduce(&local_sum, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  end_time=MPI_Wtime(); 

  MPI_Finalize() ; // let MPI finish up /
  if (myid==master)
    {
  for (i=input-resto; i < input; i++)
      sum+=array[i];
  printf ( "\n # total sum: %llu \n", sum);
    }
}
