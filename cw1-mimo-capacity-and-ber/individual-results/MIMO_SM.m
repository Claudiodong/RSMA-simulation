clc
clear
close all;
addpath("function\");
%% Parameters
Nt_MIMO = 2; % MIMO Num of transmitter antenna
Nr_MIMO = 2; % MIMO Num of receiver antenna
SNR_dB = 0:20;          % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 1e2;          % NUmber of bit generates
MAXRUN = 1e2;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;   % Number of symbol generated

%% define
Tx_symbol = zeros(Nt_MIMO,Num_symb/Nt_MIMO);
BER = zeros(1,length(SNR_dB));
BER_ZF = zeros(1,length(SNR_dB));
BER_ZF_SIC = zeros(1,length(SNR_dB));

%% QPSK modulation
% Gnerate the random bit
binary_bit = rand(1,Num_bit) > 0.5;
%QPSK 00 01 10 11
constellation = [1+1i,-1+1i,1-1i,-1-1i]/sqrt(2);

% QPSK, bit Mapping to symbol
QPSK_symbol = zeros(1,Num_bit/2);
for i = 1:Num_symb
    index = 2*binary_bit(2*i-1)+binary_bit(2*i);
    QPSK_symbol(i) = constellation(index+1);
end

% QPSK for ML decoding since transmit 2 symbols together with QPSK, then is
% have combination with 16 possible.
QPSK_ML = fML_decoding_constellation(constellation);


%% Start Monte Carlo Simulations
for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    error_sum_ML = 0;
    error_sum_ZF = 0;
    error_sum_ZF_SIC = 0;
    while (i_runs <MAXRUN)
        i_runs = i_runs+1; % count loop
        H = (randn(Nr_MIMO,Nt_MIMO)+1i*randn(Nr_MIMO,Nt_MIMO))/sqrt(2);
        noise = (randn(Nr_MIMO,Num_symb/2)+1i*randn(Nr_MIMO,Num_symb/2))/sqrt(2);
   
        %% ML deocding
        error = fMIMO_SM_with_ML(SNR(i_snr),QPSK_symbol,QPSK_ML,H,noise);
        error_sum_ML = error_sum_ML+error;

        %% Zero forcing
        error_ZF = fMIMO_SM_ZF(SNR(i_snr),QPSK_symbol,QPSK_ML,H,noise);
        error_sum_ZF = error_sum_ZF+error_ZF;

        %% nonordered ZF with SIC
        error_ZF_SIC = fMIMO_ZF_SIC_nonorder(SNR(i_snr),QPSK_symbol,constellation,H,noise);
        error_sum_ZF_SIC = error_sum_ZF_SIC + error_ZF_SIC;
    end
    % Bit error rate (BER) calculation
    BER(i_snr) = error_sum_ML/MAXRUN/Num_bit;
    BER_ZF(i_snr) = error_sum_ZF/MAXRUN/Num_bit;  
    BER_ZF_SIC(i_snr) = error_sum_ZF_SIC/MAXRUN/Num_bit;  
end
%%
figure() 
semilogy(SNR_dB,BER,"k--x","LineWidth",1);
hold on
semilogy(SNR_dB,BER_ZF,"b--x","LineWidth",1);
hold on
semilogy(SNR_dB,BER_ZF_SIC,"m--x","LineWidth",1);
legend("MIMO ML","MIMO ZF","MIMO ZF SIC nonorder")
xlabel("SNR in dB");ylabel("BER")
grid on;
% array gain arrow
annotation("textarrow",[0.5,0.15],[0.75,0.75],"String","Array Gain (Shift in SNR)");
% diversity gain arrow
annotation("textarrow",[0.83,0.7],[0.6,0.2],"String","Diversity Gain (Increase in Gradient)");
xlabel("SNR in dB");ylabel("BER");
title("2X2 MIMO System with Spatial Multiplexing and ML,ZF,SIC receiver")
%%
BER_ML = BER;

save("SM_data.mat","BER_ML","BER_ZF","BER_ZF_SIC")
