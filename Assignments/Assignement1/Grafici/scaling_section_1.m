T_comp=2*10^-9;
T_read=10^-4;
T_comm=10^-6;
np=1000;
Time_1=zeros(1, np);
Speedup_1=zeros(1, np);
N=1000000;
    for p=1:np
        Time_1(p)=T_comp*(p-1+N/p)+p*T_read+(p-1)*T_comm;
        Speedup_1(p)=Time_1(1)/Time_1(p);
    end
plot(Speedup_1)
hold on


Time_2=zeros(1, np);
Speedup_2=zeros(1, np);
N=10000000;
    for p=1:np
        Time_2(p)=T_comp*(p-1+N/p)+p*T_read+(p-1)*T_comm;
        Speedup_2(p)=Time_2(1)/Time_2(p);
    end
plot(Speedup_2, 'g')
hold on

Time_3=zeros(1, np);
Speedup_3=zeros(1, np);
N=100000000;
    for p=1:np
        Time_3(p)=T_comp*(p-1+N/p)+p*T_read+(p-1)*T_comm;
        Speedup_3(p)=Time_3(1)/Time_3(p);
    end
plot(Speedup_3, 'k')
hold on

Time_4=zeros(1, np);
Speedup_4=zeros(1, np);
N=1000000000;
    for p=1:np
        Time_4(p)=T_comp*(p-1+N/p)+p*T_read+(p-1)*T_comm;
        Speedup_4(p)=Time_4(1)/Time_4(p);
    end
plot(Speedup_4, 'r')
hold on

grid on
ylim([0,71])
title('Strong Scaling')
xlabel('Number of Processors')
ylabel('Speedup')
legend({'N=1.000.000','N=10.000.000', 'N=100.000.000', 'N=1.000.000.000'},'Location','northeast')