#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>

//binary search function

int mybsearch(int *data, int M, int N, int Key)
{
  int register low = M;
  int register high = N;
  int register mid;

  mid = (low + high) / 2;
  while(low <= high) {     

     
    if(data[mid] < Key)
      low = mid + 1; 
    else if(data[mid] > Key)
      high = mid - 1;
    else 
      return mid;

    mid = (low + high) / 2;
  }
  return -1;
}


#define N_DEFAULT  (100000000)
#define N_search_DEFAULT (N_DEFAULT / 5)


int main(int argc, char **argv)
{
  int N, Nsearch, i;
  int *data, *search;
  
  int r, tot, nthreads;
  
  int found = 0;

  //check for some inputs, otherwise use default variables
  if(argc > 1)
    N = atoi( *(argv+1) );
  else
    N = N_DEFAULT;

  if(argc > 2)
    Nsearch = atoi ( *(argv + 2) );
  else
    Nsearch = N_search_DEFAULT;

  //if openMP is not defined use the serial algorithm
  
#if !defined(_OPENMP)

  printf("performing %d lookups on %d data..\n", Nsearch, N);

  printf("set-up data.."); fflush(stdout);
  data = malloc(N * sizeof(int));
  for (i = 0; i < N; i++)
    data[i] = i;

  printf(" set-up lookups.. "); fflush(stdout);  
  search = malloc(Nsearch * sizeof(int));
  srand(time(NULL));
  for (i = 0; i < Nsearch; i++)
    search[i] = rand() % N;

  printf("\nstart cycle.. "); fflush(stdout);
    
  for (i = 0; i < Nsearch; i++)
    if( mybsearch(data, 0, N, search[i]) >= 0)
      found++;
  
#else

//pass the pointer as firstprivate to make each thread allocate the memory it needs
#pragma omp parallel firstprivate(data)  
  {

    //have a single thread allocate and initialise the array of keys
#pragma omp single
    {
      nthreads = omp_get_num_threads();

      //make sure that the problem size is divisible by the number of threads
      r=N%nthreads;
      tot=N+r;
      
      search = malloc(Nsearch * sizeof(int));
      srand(time(NULL));
      for (i = 0; i < Nsearch; i++)
	search[i] = rand() % N;
      printf("performing %d lookups on %d data..\n", Nsearch, N);
      printf("number of threads %d \n", nthreads);
    }
    
    int my_thread_id = omp_get_thread_num();

    //each thread calculates where to allocate the memory in order to avoid clashes but still have the array contiguous in memory
    //Then each thread allocates and inizialises its own memory
    data+=my_thread_id*tot/nthreads;
    data = malloc((tot/nthreads) * sizeof(int));
    
#pragma omp parallel for
    //notice that this way there are some extra elements which are inizialised,
    //but since the keys are all less than N by construction this does not affect the correctness of the code
    //If that should change one must modify the code accordingly
    //for example by making the last thread initialise the extra elements to 0 or -1 or whatever
    for (i = 0; i < tot/nthreads; i++)
      data[i] = tot*my_thread_id/nthreads+i;


    //each thread searches for the key in its own dataset
#pragma omp parallel for schedule(dynamic)    
    for (i = 0; i < Nsearch; i++)
      if( mybsearch(data, 0, tot/nthreads, search[i]) >= 0)	      
#pragma omp atomic update
	found++;  
    free(data);
  }
#endif
  
  printf("found %d values out of %d in the given dataset \n", found, Nsearch); 
  free(search);
  return 0;
}
