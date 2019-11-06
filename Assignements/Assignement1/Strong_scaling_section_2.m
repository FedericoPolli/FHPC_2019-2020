%10000000 di iterazioni
serial_time=0.25;
elapsed_time_n1=[1.62, 1.59, 1.68, 1.85, 1.89];
walltimes_n1=[0.1, 0.05, 0.028, 0.015, 0.015];
num_procs=[2, 4, 8, 16, 20];
walltime_speedup_n1=serial_time./walltimes_n1;
elapsed_speedup_n1=serial_time./elapsed_time_n1;
scatter(num_procs, elapsed_speedup_n1, 'r', 'filled')
hold on
scatter(num_procs, walltime_speedup_n1, 'b', 'filled')




ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on
title('Strong scaling of MPI pi.c for N=10.000.000')
xlabel('Number of Processors')
ylabel('Speedup')
legend({'Elapsed Time', 'WallTime'}, 'Location','northwest')