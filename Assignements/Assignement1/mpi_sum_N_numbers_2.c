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
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  FILE* in_file = fopen("input_parallel_sum.txt", "r"); 
         
  if (! in_file )
    {  
      printf("Error reading file\n"); 
      exit(-1); 
    } 

  fscanf(in_file, "%llu", &input );

  
  unsigned int* array=(unsigned  int*) calloc(input, sizeof(unsigned int));
  if (array == NULL)
    {
      printf("Error allocating memory");
      MPI_Finalize();
    }

  unsigned int i;
  unsigned int N = input/numprocs;
  unsigned int resto = input%numprocs;
  
  for (i = 0; i < N; i++)
    array[i] =myid*N + i + 1;
  
  local_sum=0;
  for (i=0; i<N ; i++) {
    local_sum+=array[i];
  }
  if (myid!=master)
    {
      printf ("\n # partial sum on processor %d: %llu\n", myid, local_sum); 
    }
  else
    {
      for (i=input-resto+1; i <= input; i++)
	local_sum+=i;
      printf("\n # partial sum on master processor: %llu\n", local_sum);
    }
  MPI_Reduce(&local_sum, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, master, MPI_COMM_WORLD);
  MPI_Finalize() ; // let MPI finish up /
  if (myid==master)
      printf ( "\n # total sum: %llu \n", sum);
}
