%10000000 di iterazioni
serial_time_n_1=0.19;
elapsed_time_n1=[0.19, 1.62, 1.59, 1.68, 1.85, 1.88];
walltimes_n1=[0.19, 0.1, 0.05, 0.028, 0.015, 0.011];
num_procs=[1, 2, 4, 8, 16, 20];
walltime_speedup_n1=serial_time_n_1./walltimes_n1;
elapsed_speedup_n1=serial_time_n_1./elapsed_time_n1;
%scatter(num_procs, elapsed_time_n1, 'r', 'filled')
hold on
%scatter(num_procs, walltimes_n1, 'b', 'filled')


%100000000 di iterazioni
serial_time_n_2=1.95;
elapsed_time_n2=[1.95, 2.54, 2.03, 1.89, 1.89, 1.97];
walltimes_n2=[1.95, 1.02, 0.51, 0.27, 0.14, 0.12];
walltime_speedup_n2=serial_time_n_2./walltimes_n2;
elapsed_speedup_n2=serial_time_n_2./elapsed_time_n2;
%scatter(num_procs, elapsed_time_n2, 'r', 'filled')
hold on
%scatter(num_procs, walltimes_n2, 'b', 'filled')


%10000000 di iterazioni
serial_time_n_3=19.54;
elapsed_time_n3=[19.54, 11.74, 6.66, 4.31, 3.20, 2.96];
walltimes_n3=[19.54, 10.2, 5.12, 2.71, 1.45, 1.15];
walltime_speedup_n3=serial_time_n_3./walltimes_n3;
elapsed_speedup_n3=serial_time_n_3./elapsed_time_n3;
scatter(num_procs, elapsed_speedup_n3, 'r', 'filled')
hold on
scatter(num_procs, walltime_speedup_n3, 'b', 'filled')



grid on
title('Strong scaling of MPI pi.c for N=1.000.000.000')
xlabel('Number of Processors')
ylabel('Speedup')
legend({'Elapsed Time', 'Walltime'}, 'Location','northeast')