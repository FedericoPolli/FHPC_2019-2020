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
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);;
  unsigned int N = input/numprocs;
  unsigned  int* array=(unsigned int*) calloc(N, sizeof(unsigned int));
  if (array == NULL)
    {
      printf("Error allocating memory");
      MPI_Finalize();
    }
  unsigned int i;
  unsigned int resto = input%numprocs;
  for (i = 0; i < N; i++)
    array[i] =myid*N + i + 1;
  // take time of processors after initial I/O operation
  start_time = MPI_Wtime();
  local_sum=0;
  for (i=0; i<N ; i++) {
    local_sum+=array[i];
  }

  if (myid != master)
    {
    MPI_Ssend(&local_sum , 1 ,MPI_LONG_LONG, master , tag ,MPI_COMM_WORLD) ;
    end_time=MPI_Wtime();
    printf ( "\n # partial sum on processor %i: %llu \n", myid, local_sum);
    printf ( "\n # walltime on processor %i : %10.8f \n\n",myid, end_time - start_time );
  }
  else
    {  
      for (i=input-resto+1; i <= input; i++)
	local_sum+=i;
      sum = local_sum;
      printf ( "\n # partial sum on master processor: %lld \n", local_sum);
      for (proc=1; proc<numprocs ; proc++) {
	MPI_Recv(&local_sum,1,MPI_LONG_LONG,proc,tag,MPI_COMM_WORLD,&status ) ;
	sum += local_sum ;
      }
      end_time=MPI_Wtime(); 
      printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ;
      printf( "\n");
      printf ( "\n # total sum: %llu \n", sum);
      printf("\n");
    }
  MPI_Finalize() ; // let MPI finish up / 
}
