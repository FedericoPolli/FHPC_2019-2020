#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <omp.h>


#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)


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

#define N_DEFAULT  (1024*1024*128)
#define N_search_DEFAULT (N_DEFAULT / 10)

int main(int argc, char **argv)
{
  int N, Nsearch, i;
  int *data, *search;
  
  int r, tot, nthreads;
  
  int found = 0;

  if(argc > 1)
    N = atoi( *(argv+1) );
  else
    N = N_DEFAULT;

  if(argc > 2)
    Nsearch = atoi ( *(argv + 2) );
  else
    Nsearch = N_search_DEFAULT;

#if !defined(_OPENMP)

  double tstart, tstop;
  struct timespec ts;
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
  tstart = CPU_TIME;
    
  for (i = 0; i < Nsearch; i++)
    if( mybsearch(data, 0, N, search[i]) >= 0)
      found++;

#else

#pragma omp parallel
  {
#pragma omp master
    {
      nthreads = omp_get_num_threads();
      r=N%nthreads;
      tot=N+r;
      search = malloc(Nsearch * sizeof(int));
      srand(time(NULL));
      data = malloc(tot * sizeof(int));
      /* for (i = 0; i < Nsearch; i++) */
      /* 	search[i] = rand() % N;  */

      /* for (i = 0; i < tot; i++) */
      /* 	data[i] = i; */
      printf("performing %d lookups on %d data..\n", Nsearch, N);
      printf("number of threads %d \n", nthreads);
    }

    #pragma omp barrier
    
    int my_thread_id = omp_get_thread_num();
    
#pragma omp parallel for
    for (i = 0; i < Nsearch; i++)
      search[i] = rand() % N;

#pragma omp parallel for
    for (i = 0; i < tot; i++)
      data[i] = i;
    
#pragma omp parallel for schedule(dynamic)
    for (i = 0; i < Nsearch; i++)
      if( mybsearch(data, my_thread_id*tot/nthreads, (my_thread_id+1)*tot/nthreads-1, search[i]) >= 0)
#pragma omp atomic
	found++;
  }
#endif
  
  printf("found: %d\n", found); 
  free(data);
  free(search);

  return 0;
}
