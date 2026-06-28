clc
clear
close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
LW = 2;

load("Nr2Q2.mat");
NT4NR2Q2 = Avg_C_sort(1,:);
NT16NR2Q2 = Avg_C_sort(2,:);
NT64NR2Q2 = Avg_C_sort(3,:);

load("Nr2Q4.mat");
NT16NR2Q4 = Avg_C_sort(1,:);
NT64NR2Q4 = Avg_C_sort(2,:);

load("Nr1Q2.mat")
NT4Q2 = Avg_C_Nt(1,:);
NT16Q2 = Avg_C_Nt(2,:);
NT64Q2 = Avg_C_Nt(3,:);
load("Nr1Q4.mat")
NT4Q4 = Avg_C_Nt(1,:);
NT16Q4 = Avg_C_Nt(2,:);
NT64Q4 = Avg_C_Nt(3,:);

CDF = (1:length(NT64NR2Q4))./length(NT16NR2Q4)*100;

figure()
subplot(2,1,1)
plot(NT4NR2Q2,CDF,'r--',"LineWidth",LW);
hold on;
plot(NT4Q2,CDF,'r-',"LineWidth",LW);
hold on;
plot(NT16NR2Q2,CDF,'b--',"LineWidth",LW);
hold on;
plot(NT16Q2,CDF,'b-',"LineWidth",LW);
hold on;
plot(NT64NR2Q2,CDF,'k--',"LineWidth",LW);
hold on;
plot(NT64Q2,CDF,'k-',"LineWidth",LW);
grid on;
legend("$N_t=4,N_r=2$","$N_t=4,N_r=1$","$N_t=16,N_r=2$","$N_t=16,N_r=1$","$N_t=64,N_r=2$","$N_t=64,N_r=1$");
xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]");
title("Average User Rate CDF under different $N_t,N_r$ and $Q=2$");

subplot(2,1,2)
plot(NT16NR2Q4,CDF,'b--',"LineWidth",LW);
hold on;
plot(NT16Q4,CDF,'b-',"LineWidth",LW);
hold on;
plot(NT64NR2Q4,CDF,'k--',"LineWidth",LW);
hold on;
plot(NT64Q4,CDF,'k-',"LineWidth",LW);
grid on;
legend("$N_t=16,N_r=2$","$N_t=16,N_r=1$","$N_t=64,N_r=2$","$N_t=64,N_r=1$");
xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]");
title("Average User Rate CDF under different $N_t,N_r$ and $Q=4$");


figure()
plot(NT4NR2Q2,CDF,"r-","LineWidth",LW)
hold on;
plot(NT16NR2Q2,CDF,"b-","LineWidth",LW)
hold on;
plot(NT64NR2Q2,CDF,"k-","LineWidth",LW)
grid on;
plot(NT16NR2Q4,CDF,"b--","LineWidth",LW)
hold on;
plot(NT64NR2Q4,CDF,"k--","LineWidth",LW)
legend("$N_t=4,Q=2$","$N_t=16,Q=2$","$N_t=64,Q=2$","$N_t=16,Q=4$","$N_t=64,Q=4$");
xlabel("Average User Rate [bps/Hz]");ylabel("Cumulative Distribution Function [\%]");
title("Average User Rate CDF with $N_r=2$ under different $N_t$ and $Q$");
xlim([0 20])

figure()
plot(NT16NR2Q4,CDF,'b--',"LineWidth",LW);
hold on;
plot(NT16Q2,CDF,'b-',"LineWidth",LW);
hold on;
plot(NT64NR2Q4,CDF,'k--',"LineWidth",LW);
hold on;
plot(NT64Q2,CDF,'k-',"LineWidth",LW);
grid on;
legend("$N_t=16,N_r=2$","$N_t=16,N_r=1$","$N_t=64,N_r=2$","$N_t=64,N_r=1$");
xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]");
title("Average User Rate CDF under different $N_t,N_r$ and $Q=4$");





