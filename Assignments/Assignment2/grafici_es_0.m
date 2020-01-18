parallel_times=[15.59, 15.33, 16.27, 17.32, 17.33, 17.34];
average_times=[14.91, 15.33, 16.27, 17.31, 17.30, 17.31];
num_threads=[1, 2, 4, 8, 16, 20];


hold on
scatter(num_threads, parallel_times, 'filled', 'b')
scatter(num_threads, average_times, 'filled', 'r')



axis([0 20 10 30])
 grid on
 xlabel('number of threads')
 ylabel('times (in seconds)')
 legend('elapsed time', "average threads' time")
 title('Weak scaling for I_{max} = 65535, with N_{pixel}/N_{thread}=500.000')
