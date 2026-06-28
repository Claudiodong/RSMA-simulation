clc
clear
close all;
addpath("function\");
%% Parameters
SNR_dB = -10:1:20;      % SNR in dB
SNR = 10.^(SNR_dB./10); % SNR in decimal
MAXRUNS = 1e4;          % The number of runs for MOnte carlo 
Nt_MIMO = [2,4,2];      % The number of transmitter antenna
Nr_MIMO = [2,2,4];      % The number of receiver antenna
Pt = 1.0;               % Transmit power 
p_tolerance = 0.00001;  % power constraint tolerance

%% predefine 
avg_power_channel_2x2_CSIT = zeros(Nt_MIMO(1),length(SNR_dB));
avg_R_2x2_CSIT = zeros(1,length(SNR_dB));
avg_R_4x2_CSIT = zeros(1,length(SNR_dB));
avg_R_2x4_CSIT = zeros(1,length(SNR_dB));
avg_R_2x2_CDIT = zeros(1,length(SNR_dB));
avg_R_4x2_CDIT = zeros(1,length(SNR_dB));
avg_R_2x4_CDIT = zeros(1,length(SNR_dB));


%% Start simulation
for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    power_sum_monte_2x2_CSIT = 0;
    % CSIT
    C_sum_2x2_CSIT = 0;
    C_sum_4x2_CSIT = 0;
    C_sum_2x4_CSIT = 0;
    % CDIT
    C_sum_2x2_CDIT = 0;
    C_sum_4x2_CDIT = 0;
    C_sum_2x4_CDIT = 0;
    % Start monte carlo 
    while(i_runs<MAXRUNS)
        i_runs = i_runs+1;
        % Generate channel

        % 2X2
        H_2x2 = (randn(Nr_MIMO(1),Nt_MIMO(1))+1i*randn(Nr_MIMO(1),Nt_MIMO(1)))/sqrt(2);
        % 4X2
        H_4x2 = (randn(Nr_MIMO(2),Nt_MIMO(2))+1i*randn(Nr_MIMO(2),Nt_MIMO(2)))/sqrt(2);
        % 2X4
        H_2x4 = (randn(Nr_MIMO(3),Nt_MIMO(3))+1i*randn(Nr_MIMO(3),Nt_MIMO(3)))/sqrt(2);


        [~,lambda_2x2,~] = svd(H_2x2);
        [~,lambda_2x4,~] = svd(H_2x4);
        [~,lambda_4x2,~] = svd(H_4x2);

        D = rank(lambda_4x2);
        B = rank(lambda_2x4);
        C = rank(lambda_2x2);             
        %% Water-Filling Algorithm
        [power_2x2_CSIT,Q_2x2_CSIT] = f_iter_water_filling(SNR(i_snr),H_2x2,Nr_MIMO(1),Pt,p_tolerance);

        [power_4x2_CSIT,Q_4x2_CSIT] = f_iter_water_filling(SNR(i_snr),H_4x2,Nr_MIMO(2),Pt,p_tolerance);

        [power_2x4_CSIT,Q_2x4_CSIT] = f_iter_water_filling(SNR(i_snr),H_2x4,Nt_MIMO(3),Pt,p_tolerance);

        % Accumulate the power
        power_sum_monte_2x2_CSIT = power_2x2_CSIT +power_sum_monte_2x2_CSIT;
        
        %% Channel Capacity

        % CSIT
        C_2x2_CSIT = real(log2(det(eye(Nr_MIMO(1))+SNR(i_snr)*H_2x2*Q_2x2_CSIT*H_2x2')));
        C_4x2_CSIT = real(log2(det(eye(Nr_MIMO(2))+SNR(i_snr)*H_4x2*Q_4x2_CSIT*H_4x2')));
        C_2x4_CSIT = real(log2(det(eye(Nr_MIMO(3))+SNR(i_snr)*H_2x4*Q_2x4_CSIT*H_2x4')));
        % Sum of all data rate
        C_sum_2x2_CSIT = C_sum_2x2_CSIT + C_2x2_CSIT;
        C_sum_4x2_CSIT = C_sum_4x2_CSIT + C_4x2_CSIT;
        C_sum_2x4_CSIT = C_sum_2x4_CSIT + C_2x4_CSIT;

        % CDIT, uniform power allocation on each channel
        C_2x2_CDIT = real(log2(det(eye(Nr_MIMO(1))+(SNR(i_snr)/Nt_MIMO(1))*(H_2x2*H_2x2'))));
        C_4x2_CDIT = real(log2(det(eye(Nr_MIMO(2))+(SNR(i_snr)/Nt_MIMO(2))*(H_4x2*H_4x2'))));
        C_2x4_CDIT = real(log2(det(eye(Nr_MIMO(3))+(SNR(i_snr)/Nt_MIMO(3))*(H_2x4*H_2x4'))));
        % Sum of all data rate
        C_sum_2x2_CDIT = C_sum_2x2_CDIT + C_2x2_CDIT;
        C_sum_4x2_CDIT = C_sum_4x2_CDIT + C_4x2_CDIT;
        C_sum_2x4_CDIT = C_sum_2x4_CDIT + C_2x4_CDIT;


    end
    % check the average power allocated to the channels for each SNR value
    avg_power_channel_2x2_CSIT(:,i_snr) = power_2x2_CSIT/MAXRUNS;
    % checking the sum of power on channels is equal to 1.
    A = sum(avg_power_channel_2x2_CSIT);
    % or we could use
    % trace(Q)

    % Average Sum Rate
    avg_R_2x2_CSIT(i_snr) = C_sum_2x2_CSIT/MAXRUNS;
    avg_R_4x2_CSIT(i_snr) = C_sum_4x2_CSIT/MAXRUNS;
    avg_R_2x4_CSIT(i_snr) = C_sum_2x4_CSIT/MAXRUNS;

    % Average Sum Rate
    avg_R_2x2_CDIT(i_snr) = C_sum_2x2_CDIT/MAXRUNS;
    avg_R_4x2_CDIT(i_snr) = C_sum_4x2_CDIT/MAXRUNS;
    avg_R_2x4_CDIT(i_snr) = C_sum_2x4_CDIT/MAXRUNS;
end

%% PLotting the graph

figure() % figure 9 in report
plot(SNR_dB,avg_R_2x2_CSIT,"b-x","LineWidth",1.5);
hold on;
plot(SNR_dB,avg_R_4x2_CSIT,"r-o","LineWidth",1.5);
hold on;
plot(SNR_dB,avg_R_2x4_CSIT,"-+","LineWidth",1.5,"Color","#4DBEEE");
hold on
plot(SNR_dB,avg_R_2x2_CDIT,"b--x","LineWidth",1.5);
hold on;
plot(SNR_dB,avg_R_4x2_CDIT,"r--o","LineWidth",1.5);
hold on;
plot(SNR_dB,avg_R_2x4_CDIT,"--+","LineWidth",1.5,"Color","#4DBEEE");
xlabel("SNR [dB]");ylabel("Ergodic Capacity [bps/Hz]");
title("The Ergodic Capacity of CSIT and CDIT in different MIMO system")
legend("2X2 CSIT","4X2 CSIT","2X4 CSIT","2X2 CDIT","4X2 CDIT","2X4 CDIT","Location","northwest")
grid on;

inserposition = [0.15,0.35,0.25,0.25];
ins = axes('Position',inserposition);
inserrange = 1:4;
hold on;
plot(ins,SNR_dB(inserrange),avg_R_2x2_CSIT(inserrange),"b-x","LineWidth",1.5);
hold on
plot(ins,SNR_dB(inserrange),avg_R_2x2_CDIT(inserrange),"b--x","LineWidth",1.5);
hold on
plot(ins,SNR_dB(inserrange),avg_R_4x2_CSIT(inserrange),"r-o","LineWidth",1.5);
hold on
plot(ins,SNR_dB(inserrange),avg_R_4x2_CDIT(inserrange),"r--o","LineWidth",1.5);
hold on
plot(ins,SNR_dB(inserrange),avg_R_2x4_CSIT(inserrange),"-+","LineWidth",1.5,"Color","#4DBEEE");
hold on
plot(ins,SNR_dB(inserrange),avg_R_2x4_CDIT(inserrange),"--+","LineWidth",1.5,"Color","#4DBEEE");
box on;
% %
% figure()
% plot(SNR_dB,avg_R_2x2_CSIT,"b-o","LineWidth",1.5);
% hold on;
% plot(SNR_dB,avg_R_2x2_CDIT,"r--o","LineWidth",1.5);
% xlabel("SNR [dB]");ylabel("Ergodic Capacity [bps/Hz]");
% title("The Ergodic Capacity of CSIT and CDIT in 2x2 MIMO system")
% legend("2X2 CSIT","2X2 CDIT")
% grid on;
% 
% 
% figure()
% plot(SNR_dB,avg_R_4x2_CSIT,"b-o","LineWidth",1.5);
% hold on;
% plot(SNR_dB,avg_R_4x2_CDIT,"r--o","LineWidth",1.5);
% xlabel("SNR [dB]");ylabel("Ergodic Capacity [bps/Hz]");
% title("The Ergodic Capacity of CSIT and CDIT in 4x2 MIMO system")
% legend("4X2 CSIT","4X2 CDIT")
% grid on;
% 
% figure()
% plot(SNR_dB,avg_R_2x4_CSIT,"b-o","LineWidth",1.5);
% hold on;
% plot(SNR_dB,avg_R_2x4_CDIT,"r--o","LineWidth",1.5);
% xlabel("SNR [dB]");ylabel("Ergodic Capacity [bps/Hz]");
% title("The Ergodic Capacity of CSIT and CDIT in 2x4 MIMO system")
% legend("2X4 CSIT","2X4 CDIT")
% grid on;
% %

