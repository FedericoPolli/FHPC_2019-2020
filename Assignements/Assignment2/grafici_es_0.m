serial_time=8.70;
%parallel_times_touch_first=[9.74, 8.35, 7.92, 7.84, 7.91, 7.94];
parallel_times_touch_all=[5.5, 3.7, 2.9, 2.7, 2.55, 2.4];
%speedup_touch_first=serial_time./parallel_times_touch_first
speedup_touch_all=serial_time./parallel_times_touch_all
%estimated_serial_part_touch_first=zeros(1, 6);
%estimated_serial_part_touch_all=zeros(1, 6);
num_threads=[2, 4, 8, 12, 16, 20];
%for i=1:6
  %  estimated_serial_part_touch_first(i)=(1/speedup_touch_first(i)-1/num_threads(i))/(1-1/num_threads(i));
    %estimated_serial_part_touch_all(i)=(1/speedup_touch_all(i)-1/num_threads(i))/(1-1/num_threads(i));
%end

%scatter([1, num_threads], [serial_time, parallel_times_touch_first], 'filled')

%scatter([1, num_threads], [serial_time, parallel_times_touch_all], 'filled', 'r')

 %scatter(num_threads, speedup_touch_first, 'filled', 'r')
% hold on
scatter(num_threads, speedup_touch_all, 'filled', 'b')


%estimated_serial_part_touch_first
%estimated_serial_part_touch_all
 axis([0 20 0 5])
 grid on
 xlabel('number of threads')
 ylabel('speedup')
 legend('binary search speedup')
