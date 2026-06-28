clc
clear
close;
addpath("function\");
%% Parameters
Nt = 2; % SISO Num of transmitter antenna
Nr = 1; % SISO Num of receiver antenna

SNR_dB = 0:20;          % SNR value in dB
SNR = 10.^(SNR_dB./10); % SNR value in decimal
Num_bit = 1e3;          % NUmber of bit generates
MAXRUN = 1e4;           % NUmber of loop for monte carlo
Num_symb = Num_bit/2;   % Number of symbol generated

%% paraeter define
BER_ala = zeros(1,length(SNR_dB));
Alamouti_QPSK_symbol = zeros(2,Num_symb);

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

% Alamouti code space time coding
% [c1,c2;-conj(c2),conj(c1)]
% loop for each Tx transmit different symbol at the same time, e.g c1 and
% c2 at Tx1 and Tx2 at 1st symbol time
for i_a = 1:Num_symb/2
    c1 = QPSK_symbol((2*i_a)-1);
    c2 = QPSK_symbol(2*i_a);
    A = (2*i_a)-1;
    B = 2*i_a;
    Alamouti_QPSK_symbol(1,A) = c1 ;
    Alamouti_QPSK_symbol(1,B) = -c2';
    Alamouti_QPSK_symbol(2,A) = c2 ;
    Alamouti_QPSK_symbol(2,B) = c1';
end

%% Start Monte Carlo Simulations
for i_snr = 1:length(SNR)
     i_max = 0;

    error_sum_MISO = 0;
    error_sum_ala = 0;

   

    while (i_max < MAXRUN)
        i_max = i_max +1;
         % channel 
        H = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2);
        % AWGN noise
        noise = (randn(Nr,Num_symb)+1i*randn(Nr,Num_symb))/sqrt(2);
    
        % Alamouti Channel 
        H_eff = [H;conj(H(2)),-conj(H(1))]; 
        % Amplitude
        p = sqrt(SNR(i_snr));
    
        % received signal, because it is transmitted with 2 symbol rate,
        % the energy per symbol is reduced by sqrt(2)
        y = p*H*(Alamouti_QPSK_symbol./sqrt(2))+noise;
       
        Z = zeros(2,Num_symb/2);
        % received signal after using match filter by using the effective
        % channel to maximise the received signal
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
              [~,index] = min(abs(reshape_symbol(i_decode)-constellation));
              estimate_symbol(i_decode) = constellation(index);
        end
        % Error calcualtion
        error = sum(real(estimate_symbol)~=real(QPSK_symbol))+sum(imag(estimate_symbol)~=imag(QPSK_symbol));
        error_sum_MISO = error_sum_MISO+error;
        error = fMISO_with_Alamouti_and_ML(SNR(i_snr),QPSK_symbol,Alamouti_QPSK_symbol,constellation,Nt,Nr);
        error_sum(i_max,i_snr) = error;
    end
    BER_ala(i_snr) = error_sum_MISO/MAXRUN/Num_bit;
end
%%
figure()
semilogy(SNR_dB,BER_ala)
