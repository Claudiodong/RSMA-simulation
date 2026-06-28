function [BER_SISO_Theory,BER_SIMO_MRC_Theory,BER_MISO_MRT_Theory,BER_MISO_Ala_Theory]=BER_theoretical(SNR_dB_range,Nr_SIMO,Nt_MISO,Nr_MISO)
SNR_dB = 0:1:SNR_dB_range;
SNR = 10.^(SNR_dB./10);
BER_SISO_Theory = zeros(1,length(SNR_dB));
BER_SIMO_MRC_Theory = zeros(1,length(SNR_dB));
BER_MISO_MRT_Theory = zeros(1,length(SNR_dB));
BER_MISO_Ala_Theory = zeros(1,length(SNR_dB));
d_min =sqrt(2); 
    for i_snr = 1:length(SNR_dB)
        % all under ML detection if not specific method is mentioned.
        % using i.i.d Rayleight fading, and channel squared is 
        % expoential distributed
        SISO = @(s) qfunc(sqrt(SNR(i_snr).*s./2).*d_min).*exp(-s);
        BER_SISO_Theory(i_snr) = integral(SISO,0,inf);
    
        SIMO_MRC = @(s) qfunc(sqrt(SNR(i_snr).*s./2).*d_min).*(1/factorial(Nr_SIMO-1)).*(s.^(Nr_SIMO-1)).*exp(-(s));
        BER_SIMO_MRC_Theory(i_snr) = integral(SIMO_MRC,0,inf);
    
        MISO_MRT = @(s) qfunc(sqrt(SNR(i_snr).*s./2).*d_min).*(1/factorial(Nr_MISO-1)).*(s.^(Nt_MISO-1)).*exp(-(s));
        BER_MISO_MRT_Theory(i_snr) = integral(MISO_MRT,0,inf); % should perfrom the same the MTC but reuqire perfect CSIT
    
        MISO_Ala = @(s) qfunc(sqrt(SNR(i_snr).*s./4).*d_min).*(1/factorial(Nr_SIMO-1)).*(s.^(Nr_SIMO-1)).*exp(-(s));
        BER_MISO_Ala_Theory(i_snr) = integral(MISO_Ala,0,inf);
    end
end