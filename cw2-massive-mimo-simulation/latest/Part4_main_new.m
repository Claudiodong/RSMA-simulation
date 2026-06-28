% clc
% clear
 close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
LW = 2;


load("Nr1Q2.mat")
NT4Q2 = Avg_C_Nt(1,:);
NT16Q2 = Avg_C_Nt(2,:);
NT64Q2 = Avg_C_Nt(3,:);
load("Nr1Q4.mat")
NT4Q4 = Avg_C_Nt(1,:);
NT16Q4 = Avg_C_Nt(2,:);
NT64Q4 = Avg_C_Nt(3,:);


CDF = (1:length(NT16Q4))./length(NT16Q4) *100;

figure()

subplot(3,1,1)
plot(NT4Q2,CDF,'r-',"LineWidth",LW);
hold on;
plot(NT4Q4,CDF,'b-.',"LineWidth",LW);
grid on;xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]")
legend("$N_t=4,Q=2$","$N_t=4,Q=4$")
title("Average User Rate CDF with $N_t=4$ and $Q=2,4$", 'Interpreter', 'latex');

subplot(3,1,2)
plot(NT16Q2,CDF,"r-","LineWidth",LW);
hold on;
plot(NT16Q4,CDF,"b-.","LineWidth",LW);
grid on;xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]")
legend("$N_t=16,Q=2$","$N_t=16,Q=4$");
title("Average User Rate CDF with $N_t=16$ and $Q=2,4$", 'Interpreter', 'latex');

subplot(3,1,3)
plot(NT64Q2,CDF,"r-","LineWidth",LW);
hold on;
plot(NT64Q4,CDF,"b-.","LineWidth",LW);
grid on;xlabel("Average User Rate [bps/Hz]");ylabel("CDF [\%]")
legend("$N_t=64,Q=2$","$N_t=64,Q=4$")
title("Average User Rate CDF with $N_t=64$ and $Q=2,4$", 'Interpreter', 'latex');

figure()
plot(NT4Q2,CDF,"r-","LineWidth",LW);
hold on;
plot(NT4Q4,CDF,"r-.","LineWidth",LW);
hold on;
plot(NT16Q2,CDF,"k-","LineWidth",LW);
hold on;
plot(NT16Q4,CDF,"k-.","LineWidth",LW);
hold on;
plot(NT64Q2,CDF,"b-","LineWidth",LW);
hold on;
plot(NT64Q4,CDF,"b-.","LineWidth",LW);
grid on;
xlabel("Average User Rate [bps/Hz]");
ylabel("Cumulative Distribution Function [\%]");
legend({'$N_t=4, Q=2$', '$N_t=4,Q=4$','$N_t=16,Q=2$','$N_t=16,Q=4$', '$N_t=64,Q=2$','$N_t=64,Q=4$'}, ...
    'Location', 'best');
title("Average User Rate CDF with different $N_t$ and $Q$", 'Interpreter', 'latex');


figure()
subplot(2,1,1)
plot(NT4Q2,CDF,"r-","LineWidth",LW);
hold on;
plot(NT16Q2,CDF,"k--","LineWidth",LW);
hold on;
plot(NT64Q2,CDF,"b-.","LineWidth",LW);
hold on;grid on;
xlabel("Average User Rate [bps/Hz]");
ylabel("CDF [\%]");
legend({'$N_t=4$', '$N_t=16$','$N_t=64$'});
title("Average User Rate CDF with different $N_t$ and $Q=2$", 'Interpreter', 'latex');

subplot(2,1,2)
plot(NT4Q4,CDF,"r-","LineWidth",LW);
hold on;
plot(NT16Q4,CDF,"k--","LineWidth",LW);
hold on;
plot(NT64Q4,CDF,"b-.","LineWidth",LW);
hold on;grid on;
xlabel("Average User Rate [bps/Hz]");
ylabel("CDF [\%]");
legend({'$N_t=4$', '$N_t=16$','$N_t=64$'});
title("Average User Rate CDF with different $N_t$ and $Q=4$", 'Interpreter', 'latex');





