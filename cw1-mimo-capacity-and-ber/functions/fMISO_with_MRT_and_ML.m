% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.14 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MISO system with Maximum Ratio Transmission at transmitter side
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

function error = fMISO_with_MRT_and_ML(SNR,Tx_QPSK_symbol,QPSK_constellation,Nt,Nr)

        Num_symb = length(Tx_QPSK_symbol);
        Num_bit = 2*Num_symb;
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
