clc
clear
close all;

%% under perfect CSIT MIMO system
SNR_dB = -10:20;
SNR = 10.^(SNR_dB./10);
Nt = 2;
Nr = 2;
MAXRUNS = 1e3;
p_tolerance = 1e-5;
Pt = 1;

avg_C1= zeros(1,length(SNR_dB));
avg_C2= zeros(1,length(SNR_dB));

for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    sum_C1 = 0;
    sum_C2 = 0;
    
    while(i_runs<MAXRUNS)
        i_runs = i_runs +1;

        H1 = [1,1;1,1];
        H2 = sqrt(2).*eye(Nt);
        % the noise power is including in the SNR value
        [power1,Q1] = f_iter_water_filling(SNR(i_snr),H1,Nr,Pt,p_tolerance);
        [power2,Q2] = f_iter_water_filling(SNR(i_snr),H2,Nr,Pt,p_tolerance);

        C1= real(log2(det(eye(Nr)+SNR(i_snr)*H1*Q1*H1')));
        C2= real(log2(det(eye(Nr)+SNR(i_snr)*H2*Q2*H2')));

        sum_C1 = sum_C1+C1;
        sum_C2 = sum_C2+C2;

    end
    avg_C1(i_snr) = sum_C1/MAXRUNS;
    avg_C2(i_snr) = sum_C2/MAXRUNS;
end


%%
figure()
plot(SNR_dB,avg_C1);
hold on
plot(SNR_dB,avg_C2);
legend("Rank 1 Channel","Rank 2 Channel")
grid on;
xlabel("SNR [dB]");ylabel("Channel Capticity [bps/Hz]");
title("The Achieveable Channel Capacity For Different Channel Condition for 2X2 MIMO system");
% Multiplexing gain arrow
annotation("textarrow",[0.8,0.8],[0.5,0.7],"String","Multiplexing Gain");