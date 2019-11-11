num_procs=[1, 2, 4, 8, 12, 16, 20];
elapsed_time_n1=[1.52, 1.53, 1.57, 1.62, 1.71, 1.82, 1.85];

scatter(num_procs, elapsed_time_n1, 'filled')

grid on
ylim([0, 2])
title('Parallel overhead of MPI pi.c')
xlabel('Number of Processors')
ylabel('Runtime')
