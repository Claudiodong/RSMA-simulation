% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.14 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MIMO system with spatial multiplexing (SM) at transmitter and having
% Zero forcing (ZF) receiver to decode the information
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  -  Tx_QPSK_symbol = modulated symbol (1x number of symbol, complex) 
%
%  - H  = channel information (2X2 for 2x2 MIMO system, complex)
%
%  - Tx_QPSK_symbol = transmit QPSK symbol information matrix
% (2xNum_symbol/2), since it apply SM at transmitter, eg if Num)symbol =
% 1000, then it have size 2x500 for 2x2 MIMO transmission system
%
%  - QPSK_ML = QPSK consetllation point for ML decoding, 
% eg QPSK for SM with
% Nt=2, then QPSK_ML have size 2X16, since different symbols are
% transmitted at the same time at different transmitter antenna, which
% result the receiver antenna will receive a combiantion of two different
% symbols information, therefore, each symbol have its own constellation,
% Then, for QPSK(4x4=16), having 16 possible combination constellation
% points.
%
%  - H = channel
%  - noise = AWGN noise
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - error_ZF = the error produced during decoding (1x1 real)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function error_ZF = fMIMO_SM_ZF(SNR,Tx_QPSK_symbol,QPSK_ML,H,noise)
        [Nr_MIMO,Nt_MIMO] = size(H);
        % Es
        p = sqrt(SNR);
        % spatial multiplexing, transmit different symbols through differet
        % antenna to imprve system diversity. with size(Nt,symbol/Nt);
        Tx_symbol = (1/sqrt(Nt_MIMO))*reshape(Tx_QPSK_symbol,Nt_MIMO,[]);
        % recived signal
        y = p*H*Tx_symbol+noise;
        %% Zero forcing
        H_ZF = pinv(H'*H)*H';% pseudo channel
        G_ZF =sqrt(Nt_MIMO/SNR)*H_ZF; % zerod forcing channel
        y_ZF = G_ZF*y; % zero forcing received signal
        estimate_symbol_ZF = zeros(Nr_MIMO,length(y_ZF));
        % ML decoding
        for i_d = 1:length(y_ZF)
            z_ZF = sum(abs(y_ZF(:,i_d)-QPSK_ML));
            [~,index_ZF] = min(z_ZF);
            estimate_symbol_ZF(:,i_d) = QPSK_ML(:,index_ZF);
        end
        % reshape the signal in to 1x[] form for decoding
        reshpae_symbol_ZF = reshape(estimate_symbol_ZF,1,[]);
        % the number of bit error occur during decoding
        error_ZF = sum(real(reshpae_symbol_ZF)~=real(Tx_QPSK_symbol))+sum(imag(reshpae_symbol_ZF)~=imag(Tx_QPSK_symbol));
end