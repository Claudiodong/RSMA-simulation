% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.14 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Generate the ML constellation for QPSK in MIMO system
% while using spatial multiplexing at transmitter 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Input:
%  - QPSK_constellation = constellation point (1x4 for QPSK complex)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - QPSK_ML = 2X16(complex) for QPSK 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function  QPSK_ML = fML_decoding_constellation(QPSK_constellation)
    % The length of constellation
    Num_QPSK = length(QPSK_constellation);
    QPSK_ML = [];
    for i = 1:Num_QPSK
        for j = 1:Num_QPSK
            QPSK = [QPSK_constellation(i);QPSK_constellation(j)];
            QPSK_ML = [QPSK_ML,QPSK];
        end
    end
end