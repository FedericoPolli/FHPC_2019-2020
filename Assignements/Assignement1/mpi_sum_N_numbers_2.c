#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define USE MPI

int main ( int argc , char *argv[ ] )
{

  long long int sum=0, local_sum ; 
   
  double start_time, end_time;   
  int myid , numprocs , proc ;
  MPI_Status status;
  MPI_Request request;
  // master process
  int master = 0;
  int tag = 123;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if ( argc <=1) {
    fprintf (stderr , " Usage : mpi -np n %s number_of_iterations \n", argv[0] ) ;
    MPI_Finalize() ;
    exit(-1) ;
  }

  unsigned long long int* array=(unsigned long long int*) calloc(atoll(argv[1]), sizeof(unsigned long long int));
  unsigned long long int i;
  long long int resto = atoll(argv[1])%numprocs;
  for (i = 0; i < atoll(argv[1]); i++)
    array[i] = i + 1;
  long long int N = atoll(argv[1])/numprocs;
// take time of processors after initial I/O operation
  start_time = MPI_Wtime();
  local_sum=0;
  for (i=0; i<N ; i++) {
    local_sum+=array[myid+i*numprocs];
  }
  MPI_Reduce(&local_sum, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  end_time=MPI_Wtime(); 

  MPI_Finalize() ; // let MPI finish up /
  for (i=atoll(argv[1])-resto; i < atoll(argv[1]); i++)
      sum+=array[i];
  if (myid==master)
  printf ( "\n # total sum: %lld \n", sum); 
}
