% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.14 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MISO system with Alamouti space time coding to improve diversity
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  - Tx_QPSK_symbol = modulated symbol (1x number of symbol, complex) 
%
%  - Alamouti_QPSK_symbol = modulated symbol for Alamouti space time coding
%  - transmission (1x number of symbol, complex) 
%
%  - QPSK_constellation = QPSK constellation point (1x4,complex)
%
%  - Nt = number of transmitter antenna (double 1x1)
%  - Nr = number of receiver antenna (double 1x1)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - error = the error produced during decoding (1x1 real)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function error = fMISO_with_Alamouti_and_ML(SNR,Tx_QPSK_symbol,Alamouti_QPSK_symbol,QPSK_constellation,Nt,Nr)

        Num_symb = length(Tx_QPSK_symbol);
        % channel 
        H = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2);
        % AWGN noise
        noise = (randn(Nr,Num_symb)+1i*randn(Nr,Num_symb))/sqrt(2);
    
        % Alamouti Channel 
        H_eff = [H;conj(H(2)),-conj(H(1))]; 
        % Amplitude
        p = sqrt(SNR);
    
        % received signal, because it is transmitted with 2 symbol rate,
        % the energy per symbol is reduced by sqrt(2)
        y = p*H*(Alamouti_QPSK_symbol./sqrt(2))+noise;
       
        Z = zeros(2,Num_symb/2);
        % received signal after using match filter
        for i = 1:Num_symb/2
            y1 = y((2*i)-1);
            y2 = y(2*i);
            Z(:,i) = H_eff'*[y1;conj(y2)];
        end

        % reshape the signal into 1xNumsymbol
        reshape_symbol = reshape(Z,1,[]);

        % ML decoding
        estimate_symbol = zeros(1,Num_symb);
        for i_decode = 1:Num_symb
              [~,index] = min(abs(reshape_symbol(i_decode)-QPSK_constellation));
              estimate_symbol(i_decode) = QPSK_constellation(index);
        end
        % Error calcualtion
        error = sum(real(estimate_symbol)~=real(Tx_QPSK_symbol))+sum(imag(estimate_symbol)~=imag(Tx_QPSK_symbol));
end