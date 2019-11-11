%N=10.000.000
serial_time_n1=1.55;
elapsed_time_n1=[1.55, 1.55, 1.57, 1.62, 1.83, 1.84];
num_procs=[1, 2, 4, 8, 16, 20];
elapsed_speedup_n1=serial_time_n1./elapsed_time_n1;
scatter(num_procs, elapsed_speedup_n1, 'y', 'filled')
hold on
%scatter(num_procs, walltimes_n1, 'b', 'filled')


%N=100.000.000
serial_time_n2=1.72;
elapsed_time_n2=[1.72, 1.62, 1.58, 1.65, 1.84, 1.85];
elapsed_speedup_n2=serial_time_n2./elapsed_time_n2;
scatter(num_procs, elapsed_speedup_n2, 'b', 'filled')
hold on
%scatter(num_procs, walltimes_n2, 'b', 'filled')


%N=1.000.000.000
serial_time_n3=3.88;
elapsed_time_n3=[3.88, 2.75, 2.17, 1.95, 1.89, 1.95];
elapsed_speedup_n3=serial_time_n3./elapsed_time_n3;
scatter(num_procs, elapsed_speedup_n3, 'r', 'filled')
hold on
%scatter(num_procs, walltimes_n3, 'b', 'filled')


grid on
title('Strong Scaling')
xlabel('Number of Processors')
ylabel('Speedup')
legend({'N=10.000.000', 'N=100.000.000', 'N=1.000.000.000'},'Location','northeast')