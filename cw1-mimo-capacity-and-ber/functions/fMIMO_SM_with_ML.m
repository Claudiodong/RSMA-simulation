% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.14 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MIMO system with spatial multiplexing (SM) at transmitter and having
% Maximum likelihood (ML) receiver to decode the information
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  - Tx_QPSK_symbol = modulated symbol (1x number of symbol, complex) 
%
%  - QPSK_ML = QPSK consetllation point for ML decoding, eg QPSK for SM with
% Nt=2, then QPSK_ML have size 2X16, since different symbols are
% transmitted at the same time at different transmitter antenna, which
% result the receiver antenna will receive a combiantion of two different
% symbols information, therefore, each symbol have its own constellation,
% Then, for QPSK(4x4=16), having 16 possible combination constellation
% points.
%
%  - H = Channel
%  - Noise = Noise
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - error = the error produced during decoding (1x1 real)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function error = fMIMO_SM_with_ML(SNR,Tx_QPSK_symbol,QPSK_ML,H,noise)
%        Num_symb = length(Tx_QPSK_symbol);
        [Nr_MIMO,Nt_MIMO] = size(H);
        
        % amplitude
        p = sqrt(SNR);
        % spatial multiplexing, transmit different symbols through differet
        % antenna to imprve system diversity. with size(Nt,symbol/Nt);
        Tx_symbol = (1/sqrt(Nt_MIMO))*reshape(Tx_QPSK_symbol,Nt_MIMO,[]);

        % recived signal
        y = p*H*Tx_symbol+noise;

        % ML decoding
        estimate_symbol = zeros(Nr_MIMO,length(y));
        for i_d = 1:length(y)
            % comparing one received signal to all 16 combination
            % constellation points
            z = sum(abs(y(:,i_d)-p*H*(1/sqrt(Nt_MIMO))*QPSK_ML));
            % find the minimum distance, which means the cloest
            % consetllation point that it will map to
            [~,index] = min(z);
            % estimated symbol after using ML receiver.
            estimate_symbol(:,i_d) = QPSK_ML(:,index);
        end

        % reshape the signal into an 1xNum_symbol to detect the error
        reshpae_symbol = reshape(estimate_symbol,1,[]);
        error = sum(real(reshpae_symbol)~=real(Tx_QPSK_symbol))+sum(imag(reshpae_symbol)~=imag(Tx_QPSK_symbol));

end
