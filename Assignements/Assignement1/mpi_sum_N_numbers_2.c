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

  //each processor reads input from file
  FILE* in_file = fopen("input_parallel_sum.txt", "r"); 
         
  if (! in_file )
    {  
      printf("Error reading file\n"); 
      exit(-1); 
    } 

  fscanf(in_file, "%llu", &input );

  unsigned int i;
  unsigned int N = input/numprocs;
  unsigned int resto = input%numprocs;
  
  start_time = MPI_Wtime();
  local_sum=0;
  //each processor computes a partial sum
  for (i = 0; i < N; i++)
    local_sum+=myid*N + i + 1;

  if (myid!=master)
    {
      printf ("\n # partial sum on processor %d: %llu\n", myid, local_sum); 
    }
  else
    {
      //master takes care of remaining numbers should numproc not divide input
      for (i=input-resto+1; i <= input; i++)
	local_sum+=i;
      printf("\n # partial sum on master processor: %llu\n", local_sum);
    }
  //computes total sum
  MPI_Reduce(&local_sum, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, master, MPI_COMM_WORLD);
  end_time=MPI_Wtime();
  if (myid==master){
     printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ;
    printf ( "\n # total sum: %llu \n", sum);
  }
  else {
    printf ( "\n # walltime on processor %i : %10.8f \n\n",myid, end_time - start_time );}
  MPI_Finalize() ; // let MPI finish up /

}
