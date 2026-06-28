clc
clear
close;
addpath("function\");
%%
Nr = 2; % Num of receiver antenna
Nt = 1; % Num of transmitter antenna

SNR_dB = 0:20;          % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 1e3;          % NUmber of bit generates
MAXRUN = 1e3;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;   % Number of symbol generated

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

for i_snr = 1:length(SNR_dB)
    i_runs = 0;
    error_sum = 0;
    while (i_runs <MAXRUN)
        i_runs = i_runs+1;
       
        error_SIMO = fSIMO_with_MRC_and_ML(SNR(i_snr),QPSK_symbol,constellation,Nt,Nr);
        error_sum =error_sum +error_SIMO;
    end
    % Bit error rate calculation
    BER(i_snr) = error_sum/MAXRUN/Num_bit;
end

dmin = sqrt(2);
BER_siso = [];
N = 2;
for i=1:size(SNR_dB,2)
    snr = 10.^(SNR_dB(i)./10);
    fun_siso = @(u) qfunc(sqrt(snr.*u./2) .* dmin) .* 1./factorial(N-1) .* u.^(N-1) .* exp(-u); 
    BER_siso = [BER_siso, integral(fun_siso, 0, Inf)];
end

figure()
semilogy(SNR_dB,BER);
grid on;
%%
save('SIMO.mat','BER')



