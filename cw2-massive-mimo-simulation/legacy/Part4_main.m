clc
clear
close all;

%%
% Q2(1,:) Nt=4
% Q2(2,:) Nt=16
% Q2(3,:) Nt=64
% But Q=2

load("Q2.mat")
Q2 = Avg_C_all;
%%
load("Q4.mat")
Q4 = Avg_C_all;

cdf = (1:length(Avg_C_all))./length(Avg_C_all);

figure()
plot(Q2(1,:),cdf,"r--",'LineWidth',1.5)
hold on;
plot(Q4(1,:),cdf,"r-.",'LineWidth',1.5)
hold on;
plot(Q2(2,:),cdf,"b--",'LineWidth',1.5)
hold on;
plot(Q4(2,:),cdf,"b-.",'LineWidth',1.5)
hold on;
plot(Q2(3,:),cdf,"k--",'LineWidth',1.5)
hold on;
plot(Q4(3,:),cdf,"k-.",'LineWidth',1.5)
grid on;
xlabel("Average User Rate [bps/Hz]");ylabel("Cumulative Distribution Function [%]");
legend("N_t=4, Q=2","N_t=4, Q=4","N_t=16, Q=2","N_t=16, Q=4","N_t=64, Q=2","N_t=64, Q=4")

figure()
plot(Q2(1,:),cdf,"r--",'LineWidth',1.5)
hold on;
plot(Q4(1,:),cdf,"r-.",'LineWidth',1.5)
legend("N_t=4, Q=2","N_t=4, Q=4");grid on
