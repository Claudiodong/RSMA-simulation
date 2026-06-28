clc
clear
close all;
addpath("function\");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SISO system with QPSK modulation, default using ML decoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_dB = 0:20;          % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 1e3;          % NUmber of bit generates
MAXRUN = 1e3;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;
BER = zeros(1,length(SNR_dB));
y_hat = zeros(1,Num_symb);
Nt = 1;
Nr = 1;

% Gnerate the random bit
binary_bit = rand(1,Num_bit) > 0.5;
%QPSK 00 01 10 11
constellation = [1+1i,-1+1i,1-1i,-1-1i]/sqrt(2);

% QPSK bit Mapping to symbol
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
        % error by ML decoding
        [error] = fSISO_with_ML(QPSK_symbol,SNR(i_snr),constellation,Nt,Nr);
        error_sum = error+error_sum;     
    end
    BER_SISO(i_snr) = error_sum/MAXRUN/Num_bit;
end

% dmin = sqrt(2);
% BER_siso = [];
% for i=1:size(SNR_dB,2)
%     snr = 10.^(SNR_dB(i)./10);
%     fun_siso = @(u) qfunc(sqrt(snr.*u./2) .* dmin) .* exp(-u); 
%     BER_siso = [BER_siso, integral(fun_siso, 0, Inf)];
% end

semilogy(SNR_dB,BER_SISO)
hold on
semilogy(SNR_dB,BER_siso)
grid on
%%
save('SISO.mat','BER_SISO')
