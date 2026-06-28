clc
clear
close all;
addpath("function\");
%% Parameters
% SISO
Nt_SISO = 1; % SISO Num of transmitter antenna
Nr_SISO = 1; % SISO Num of receiver antenna
% SIMO
Nt_SIMO = 1; % SIMO Num of transmitter antenna
Nr_SIMO = 2; % SIMO Num of receiver antenna
% MISO
Nt_MISO = 2; % SIMO Num of transmitter antenna
Nr_MISO = 1; % SIMO Num of receiver antenna
% MIMO
Nt_MIMO = 2; % MIMO Num of transmitter antenna
Nr_MIMO = 2; % MIMO Num of receiver antenna
% Others
SNR_dB = 0:1:20;        % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 5e3;          % NUmber of bit generates
MAXRUN = 1e5;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;   % Number of symbol generated

%%  Predefine
BER_SISO = zeros(1,length(SNR_dB));
BER_SIMO = zeros(1,length(SNR_dB));
BER_MISO = zeros(1,length(SNR_dB));
BER_ala = zeros(1,length(SNR_dB));
BER_MIMO_SM_ML = zeros(1,length(SNR_dB));
BER_MIMO_SM_ZF = zeros(1,length(SNR_dB));
BER_MIMO_SM_ZF_SIC = zeros(1,length(SNR_dB));

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

figure() % figure 1 in report
plot(real(constellation),imag(constellation),"O","MarkerSize",20);
str = {'00','01','10','11'};
text(real(constellation)-0.02,imag(constellation),str);
title("QPSK Constellation");xlabel("Real");ylabel("Imag");grid on
%% Alamouti code space time coding

% [c1,c2;-conj(c2),conj(c1)]
Alamouti_QPSK_symbol = zeros(2,Num_symb/2);
% loop for each Tx transmit different symbol at the same time, e.g c1 and
% c2 at Tx1 and Tx2 at 1st symbol time
for i_a = 1:Num_symb/2
    Alamouti_QPSK_symbol(1,i_a) = QPSK_symbol((2*i_a)-1);
    Alamouti_QPSK_symbol(2,i_a) = QPSK_symbol(2*i_a);
end

%% MIMO system perform spatial multiplexing to improve the diveristy
% QPSK for ML decoding since transmit 2 symbols together with QPSK, then is
% have combination with 16 possible.
QPSK_ML = fML_decoding_constellation(constellation);

%% Start Monte Carlo Simulations
for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    error_sum_SISO = 0;
    error_sum_SIMO = 0;
    error_sum_MISO = 0;
    error_sum_ala = 0;
    error_sum_MIMO_SM_ML = 0;
    error_sum_MIMO_SM_ZF = 0;
    error_sum_ZF_SIC = 0;
    while (i_runs <MAXRUN)
        i_runs = i_runs+1; % count loop

        %% SISO with ML
        error_SISO = fSISO_with_ML(QPSK_symbol,SNR(i_snr),constellation,Nt_SISO,Nr_SISO);
        error_sum_SISO = error_sum_SISO+error_SISO;
      
        %% SIMO with MRC with ML
        error_SIMO = fSIMO_with_MRC_and_ML(SNR(i_snr),QPSK_symbol,constellation,Nt_SIMO,Nr_SIMO);
        error_sum_SIMO =error_sum_SIMO +error_SIMO;

        %% MISO with MRT with ML
        error_MISO = fMISO_with_MRT_and_ML(SNR(i_snr),QPSK_symbol,constellation,Nt_MISO,Nr_MISO);
        error_sum_MISO =error_sum_MISO +error_MISO;

        %% MISO with Alamouti space time coding with ML
        error_Ala = fMISO_with_Alamouti_and_ML(SNR(i_snr),QPSK_symbol,Alamouti_QPSK_symbol,constellation,Nt_MISO,Nr_MISO);
        error_sum_ala = error_sum_ala+error_Ala;

        %% MIMO channel 
        H_MIMO = (randn(Nr_MIMO,Nt_MIMO)+1i*randn(Nr_MIMO,Nt_MIMO))/sqrt(2);
        % MIMO noise
        noise_MIMO = (randn(Nr_MIMO,Num_symb/2)+1i*randn(Nr_MIMO,Num_symb/2))/sqrt(2);

        %% MIMO with SM and ML
        error_MIMO_SM_ML = fMIMO_SM_with_ML(SNR(i_snr),QPSK_symbol,QPSK_ML,H_MIMO,noise_MIMO);
        error_sum_MIMO_SM_ML = error_sum_MIMO_SM_ML + error_MIMO_SM_ML;

        %% MIMO with SM and ZF
        error_MIMO_SM_ZF = fMIMO_SM_ZF(SNR(i_snr),QPSK_symbol,QPSK_ML,H_MIMO,noise_MIMO);
        error_sum_MIMO_SM_ZF = error_sum_MIMO_SM_ZF+error_MIMO_SM_ZF;

        %% MIMO with SM and ZF and SIC
        error_ZF_SIC = fMIMO_ZF_SIC_nonorder(SNR(i_snr),QPSK_symbol,constellation,H_MIMO,noise_MIMO);
        error_sum_ZF_SIC = error_sum_ZF_SIC + error_ZF_SIC;
    end
    % Bit error rate (BER) calculation
    BER_SISO(i_snr) = error_sum_SISO/MAXRUN/Num_bit;
    BER_SIMO(i_snr) = error_sum_SIMO/MAXRUN/Num_bit;
    BER_MISO(i_snr) = error_sum_MISO/MAXRUN/Num_bit;
    BER_ala(i_snr) = error_sum_ala/MAXRUN/Num_bit;
    BER_MIMO_SM_ML(i_snr) = error_sum_MIMO_SM_ML/MAXRUN/Num_bit;
    BER_MIMO_SM_ZF(i_snr) = error_sum_MIMO_SM_ZF/MAXRUN/Num_bit;
    BER_MIMO_SM_ZF_SIC(i_snr) = error_sum_ZF_SIC/MAXRUN/Num_bit;
end
%% Theoretical result for comparsion
[BER_SISO_Theory,BER_SIMO_MRC_Theory,BER_MISO_MRT_Theory,BER_MISO_Ala_Theory]=BER_theoretical(20,Nr_SIMO,Nt_MISO,Nr_MISO);
%% plotting 
figure() % figure 8 in report
semilogy(SNR_dB,BER_SISO,"k--x","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_SIMO,"r--o","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_MISO,"b--x","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_ala,"g--x","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_MIMO_SM_ML,"--x","LineWidth",1.5,'Color',[0.8500 0.3250 0.0980],"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_MIMO_SM_ZF,"--x","LineWidth",1.5,'Color',[0.4660 0.6740 0.1880],"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_MIMO_SM_ZF_SIC,"--x","LineWidth",1.5,'Color',"#4DBEEE","MarkerSize",8);
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("SISO","SIMO with MRC","MISO with MRT","MISO with Alamouti" ...
    ,"MIMO with SM and ML","MIMO with SM and ZF","MIMO with SM and nonorder-ZF-SIC","Location","southwest");
title("BER against SNR with different systems and reception stragies");
% array gain arrow
annotation("textarrow",[0.6,0.2],[0.75,0.75],"String","Array Gain");
% diversity gain arrow
annotation("textarrow",[0.83,0.7],[0.7,0.2],"String","Diversity Gain");


%%
figure()  % figure 7 in report
semilogy(SNR_dB,BER_MIMO_SM_ML,"k--x","LineWidth",1);
hold on
semilogy(SNR_dB,BER_MIMO_SM_ZF,"b--x","LineWidth",1);
hold on
semilogy(SNR_dB,BER_MIMO_SM_ZF_SIC,"m--x","LineWidth",1);
legend("MIMO ML","MIMO ZF","MIMO ZF SIC nonorder")
xlabel("SNR in dB");ylabel("BER")
grid on;
% array gain arrow
annotation("textarrow",[0.45,0.2],[0.72,0.72],"String","Array Gain (Shift in SNR)");
% diversity gain arrow
annotation("textarrow",[0.80,0.7],[0.6,0.2],"String","Diversity Gain (increase in gradient)");
xlabel("SNR in dB");ylabel("BER");
title("2X2 MIMO System with Spatial Multiplexing and ML,ZF,SIC receiver")

%% plotting 
figure() % figure 6 in report
semilogy(SNR_dB,BER_SISO,"k--x","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_SIMO,"r--o","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_MISO,"b--x","LineWidth",1.5,"MarkerSize",8);
hold on
semilogy(SNR_dB,BER_ala,"g--x","LineWidth",1.5,"MarkerSize",8);
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("SISO","SIMO with MRC","MISO with MRT","MISO with Alamouti");
title("BER against SNR with different systems");
% array gain arrow
annotation("textarrow",[0.6,0.2],[0.7,0.7],"String","Array Gain");
% diversity gain arrow
annotation("textarrow",[0.83,0.7],[0.7,0.2],"String","Diversity Gain");


%% individual comparsion for each BER in practical and theoretical

figure() % figure 2 in report
% For SISO
semilogy(SNR_dB,BER_SISO,"r-x");
hold on;
semilogy(SNR_dB,BER_SISO_Theory,"k--o");
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("Practical SISO","Theoretical SISO");
title("Comparsion of Practical and Theoretical SISO system");

%%
figure()  % figure 3 in report
% For SIMO
semilogy(SNR_dB,BER_SIMO,"b-x");
hold on;
semilogy(SNR_dB,BER_SIMO_MRC_Theory,"k--o");
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("Practical SIMO MRC","Theoretical SIMO MRC");
title("Comparsion of Practical and Theoretical SIMO MRC system");

%%
figure() % figure 4 in report
% For MISO
semilogy(SNR_dB,BER_MISO,"b-x");
hold on;
semilogy(SNR_dB,BER_MISO_MRT_Theory,"k--o");
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("Practical MISO MRT","Theoretical MISO MRT");
title("Comparsion of Practical and Theoretical MISO MRT system");

%%
figure() % figure 5 in report
% For MISO alamouti
semilogy(SNR_dB,BER_ala,"b-x");
hold on;
semilogy(SNR_dB,BER_MISO_Ala_Theory,"k--o");
grid on;xlabel("SNR [dB]");ylabel("Bit Error Rate (BER) [bps/Hz]");
legend("Practical MISO Alamouti","Theoretical MISO Alamouti");
title("Comparsion of Practical and Theoretical MISO Alamouti system");

%%
% figure()
% % SISO
% subplot(2,2,1);
% semilogy(SNR_dB,BER_SISO,"b-x");hold on;semilogy(SNR_dB,BER_SISO_Theory,"k--o");
% grid on;xlabel("SNR [dB]");ylabel("BER [bps/Hz]");title("SISO System");
% legend("Practical","Theoretical");
% % SIMO MRC
% subplot(2,2,2);
% semilogy(SNR_dB,BER_SIMO,"b-x");hold on;semilogy(SNR_dB,BER_SIMO_MRC_Theory,"k--o");
% grid on;xlabel("SNR [dB]");ylabel("BER [bps/Hz]");title("SIMO System");
% legend("Practical","Theoretical");
% % MISO MRT
% subplot(2,2,3);
% semilogy(SNR_dB,BER_MISO,"b-x");hold on;semilogy(SNR_dB,BER_MISO_MRT_Theory,"k--o");
% grid on;xlabel("SNR [dB]");ylabel("BER [bps/Hz]");title("MISO MRT System");
% legend("Practical","Theoretical");
% % MISO Alamouti
% subplot(2,2,4);
% semilogy(SNR_dB,BER_ala,"b-x");hold on;semilogy(SNR_dB,BER_MISO_Ala_Theory,"k--o");
% grid on;xlabel("SNR [dB]");ylabel("BER [bps/Hz]");title("MISO Alamouti System");
% legend("Practical","Theoretical");

