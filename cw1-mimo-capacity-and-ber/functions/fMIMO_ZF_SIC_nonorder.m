% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.17 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% MIMO system with spatial multiplexing (SM) at transmitter and having
% Zero forcing (ZF) and Successive interference canceler (SIC) to decode the information
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  - Tx_QPSK_symbol = modulated symbol (1x number of symbol, complex) 
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
%  - H_MIMO = MIMO channel (NrxNt complex)
%  - noise_MIMO = MIMO noise (NrxNum_symbol/2, complex)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - error_ZF_SIC = the error produced during decoding (1x1 real)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function error_ZF_SIC = fMIMO_ZF_SIC_nonorder(SNR,Tx_QPSK_symbol,QPSK_constellation,H_MIMO,noise_MIMO)
        [Nr_MIMO,Nt_MIMO] = size(H_MIMO);
        % number of transmit symbol without SM
        Num_symb = length(Tx_QPSK_symbol);
        % Es
        d = sqrt(SNR);
        % transmit symbol with spatial multiplexing
        Tx_symbol = (1/sqrt(Nt_MIMO))*reshape(Tx_QPSK_symbol,Nt_MIMO,[]);
        % received signal
        y = d*H_MIMO*Tx_symbol+noise_MIMO;

        % G_ZF
        H_ZF = pinv(H_MIMO'*H_MIMO)*H_MIMO';% pseudo channel
        G_ZF =sqrt(Nt_MIMO/SNR)*H_ZF; % zerod forcing channel

        % received signal with ZF
        y_ZF = G_ZF*y;

        % pre-define the SIC channel parameter
        estimate_symbol_ZF = zeros(Nr_MIMO,Num_symb/2);
        re_encode = zeros(Nt_MIMO-1,Num_symb/2);
        H_SIC = H_MIMO;
        %  decoding
        for i_Nt = 1:Nt_MIMO
            for i = 1:Num_symb/2
                % measure the distance to each constellation point
                distance = abs(y_ZF(i_Nt,i)-QPSK_constellation);
                % find the closet point and map to it
                [~,index_zf] = min(distance);
                estimate_symbol_ZF(i_Nt,i) = QPSK_constellation(index_zf);
            end
            % perfrom SIC until the last antenna
            
            for i_sic = i_Nt+1:Nt_MIMO
                % re-encode the estimated symbol, but if it is perfect
                % should use estimate_symbol_ZF or Tx_symbol
                re_encode = sqrt(SNR/Nt_MIMO)*H_MIMO(:,i_Nt)*estimate_symbol_ZF(i_Nt,:)+re_encode;
                % Perform SIC, assumed with perfect SIC
                y_SIC = y - re_encode;
                % reproduce the ZF matirx based on new channel since the first
                % symbol is substarcted and its channel.
                % Using pinv to prevent singular matrix produced
                H_SIC(:,i_Nt) = zeros(Nt_MIMO,1);
                % pesudo inverse Channel 
                H_SIC_pse = pinv(H_SIC'*H_SIC)*H_SIC';
                G_SIC =  sqrt(Nt_MIMO/SNR)*H_SIC_pse;
                % Using the new G_ZF to apply on the recived signal with SIC
                y_ZF = G_SIC*y_SIC;
            end
        end
       
        % Reshape the signal in to 1x[] form for decoding
        reshpae_symbol_ZF = reshape(estimate_symbol_ZF,1,[]);
        % The number of bit error occur during decoding
        error_ZF_SIC = sum(real(reshpae_symbol_ZF)~=real(Tx_QPSK_symbol))+sum(imag(reshpae_symbol_ZF)~=imag(Tx_QPSK_symbol));

end