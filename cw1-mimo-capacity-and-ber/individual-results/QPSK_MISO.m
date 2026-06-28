clc
clear
close;
addpath("function\");
%% Parameters
Nt_MISO = 2; % SISO Num of transmitter antenna
Nr_MISO = 1; % SISO Num of receiver antenna

SNR_dB = 0:20;          % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 1e3;          % NUmber of bit generates
MAXRUN = 1e3;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;   % Number of symbol generated

%% paraeter define
BER_MISO = zeros(1,length(SNR_dB));

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

%% Start Monte Carlo Simulations
for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    error_sum_MISO = 0;
    while (i_runs <MAXRUN)
        i_runs = i_runs+1; % count loop
        error = MISO_with_MRT_and_ML(SNR(i_snr),QPSK_symbol,constellation,Nt_MISO,Nr_MISO,Num_bit,Num_symb);
        error_sum_MISO = error_sum_MISO +error;
    end
    % Bit error rate (BER) calculation
    BER_MISO(i_snr) = error_sum_MISO/MAXRUN/Num_bit;
    
end

figure()
semilogy(SNR_dB,BER_MISO)


function error = MISO_with_MRT_and_ML(SNR,Tx_QPSK_symbol,QPSK_constellation,Nt,Nr,Num_bit,Num_symb)

        % channel
        H = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2); 
        % Noise
        noise = (randn(Nr,Num_bit/2)+1i*randn(Nr,Num_bit/2))/sqrt(2);
        % Amplitude
        p = sqrt(SNR);

        % precoder MRT design
        w = H'/norm(H);

        % received signal
        y = (p*H*w)*Tx_QPSK_symbol+noise;

        % ML decoding
         % ML decoding
        estimate_symbol = zeros(1,Num_symb);
        for i = 1:Num_symb
            % Because Z have factor out the H since the g is in terms of H
            % and the SNR is maximized.
            [~,index] = min(abs(y(i)-QPSK_constellation));
            estimate_symbol(i) = QPSK_constellation(index);
        end

         % bit error calculation in terms of real and imag
        error = sum(real(estimate_symbol)~=real(Tx_QPSK_symbol))+sum(imag(estimate_symbol)~=imag(Tx_QPSK_symbol));
end

