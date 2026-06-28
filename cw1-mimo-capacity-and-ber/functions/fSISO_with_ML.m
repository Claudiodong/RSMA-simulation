% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.22 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% SISO system with Maximum Ratio Transmission at transmitter side
% which is precoder design for the system to maximise system SNR
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  - Tx_QPSK_symbol = modulated symbol (1x number of symbol, complex) 
%
%  - QPSK_constellation = QPSK constellation point (1x4,complex)
%
%  -  Nt = number of transmitter antenna (double 1x1)
%  -  Nr = number of receiver antenna (double 1x1)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - error = the error produced during decoding (1x1 real)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [error] = fSISO_with_ML(Tx_QPSK_symbol,SNR,QPSK_constellation,Nt,Nr)
        Num_symb = length(Tx_QPSK_symbol);
        Num_bit = 2*Num_symb;
        % channel
        H = (randn(Nt,Nr)+1i*randn(Nt,Nr))/sqrt(2); 
        % Noise
        noise = (randn(Nr,Num_bit/2)+1i*randn(Nr,Num_bit/2))/sqrt(2);
        p = sqrt(SNR);
        % received signal
        y = p*H*Tx_QPSK_symbol+noise;
        % ML decoding
        y_hat =zeros(1,Num_symb);
        for i = 1:length(y)
            [~,index] = min(abs(y(i) - p*H*QPSK_constellation));   
            y_hat(i) = QPSK_constellation(index);
        end
        error = sum(real(y_hat)~=real(Tx_QPSK_symbol))+sum(imag(y_hat)~=imag(Tx_QPSK_symbol));
end