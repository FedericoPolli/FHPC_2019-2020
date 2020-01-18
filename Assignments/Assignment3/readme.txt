The file exercise1.cpp implements a code that calculates, in parallel with OpenMP, the Mandelbrot set.
It can be compiled with the -DBALANCED flag to provide a more balanced, albeit slightly more inefficient, implementation.
On Ulysses, before running the executable, the environmental variable OMP_PROC_BIND should be set to true to avoid thread migration.

The second file implemets a hybrid MPI+OpenMP version. Such code is meant to run with the number of processes that divides the input size n_x*n_y, and that also divides the total number of cores available. this way each process will spawn exactly (total number of processors)/(number of processes started) threads, so as not to waste any computational power.
Moreover the executable should be run with the flag --bind-to none, otherwise OpenMP will not work properly.

