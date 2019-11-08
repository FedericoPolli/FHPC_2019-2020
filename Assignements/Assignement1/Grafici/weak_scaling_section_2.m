%N/p=500000000
serial_time_n_1=0.98;
elapsed_time_n1=[0.98, 2.52, 2.56, 2.71, 2.85, 2.89, 2.94];
walltimes_n1=[0.98, 1.02, 1.02, 1.08, 1.15, 1.15, 1.16];
num_procs=[1, 2, 4, 8, 12, 16, 20];
walltime_speedup_n1=serial_time_n_1./walltimes_n1;
elapsed_speedup_n1=serial_time_n_1./elapsed_time_n1;
scatter(num_procs, elapsed_speedup_n1, 'r', 'filled')
hold on
scatter(num_procs, walltime_speedup_n1, 'b', 'filled')

grid on
title('Weak scaling of MPI Pi with N/nproc=500.000.000')
xlabel('Number of Processors')
ylabel('Speedup')
legend({'Elapsed Time', 'Walltime'}, 'Location','northeast')