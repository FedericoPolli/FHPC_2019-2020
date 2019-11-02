#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define USE MPI

int main ( int argc , char *argv[ ] )
{

  // number of points inside the circle
  long long int sum=0, local_sum ; 
   
  // times 
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

  if (myid == master) { //if I am the master process gather results from others
    for (i=atoll(argv[1])-resto; i < atoll(argv[1]); i++)
      local_sum+=array[i];
    sum = local_sum ;
    printf ( "\n # partial sum on master processor: %lld \n", local_sum);
    for (proc=1; proc<numprocs ; proc++) {
      MPI_Recv(&local_sum,1,MPI_LONG_LONG,proc,tag,MPI_COMM_WORLD,&status ) ;

      sum += local_sum ;
    }
    end_time=MPI_Wtime(); 
    printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ;

  }
  else {   // for all the slave processes send results to the master /

//    printf ( " Processor %d sending results = %llu to master process \n", myid, local_M) ;
//    int time_to_sleep=1*myid;
//    sleep(time_to_sleep);

    MPI_Ssend(&local_sum , 1 ,MPI_LONG_LONG, master , tag ,MPI_COMM_WORLD) ;
    end_time=MPI_Wtime();
    printf ( "\n # partial sum on processor %i: %lld \n", myid, local_sum);
    printf ( "\n # walltime on processor %i : %10.8f \n",myid, end_time - start_time ) ;
  }
  MPI_Finalize() ; // let MPI finish up /
  printf ( "\n # total sum: %lld \n", sum); 
}
