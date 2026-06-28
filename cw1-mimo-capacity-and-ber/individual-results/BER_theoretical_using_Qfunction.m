clc
%clear
close all;

SNR_dB = 0:2:20;
SNR = 10.^(SNR_dB./10);
d_min =sqrt(2); % because it is QPSK, then the minimum distance between the
% normalised points is sqrt(2);

Nr_SIMO = 2;
Nt_MISO = 2;
Nr_MIMO = 2;
Nt_MIMO = 2;

%QPSK 00 01 10 11
constellation = [1+1i,-1+1i,1-1i,-1-1i]/sqrt(2);


BER_SISO_Theory = zeros(1,length(SNR_dB));
BER_SIMO_MRC_Theory = zeros(1,length(SNR_dB));
BER_MISO_MRT_Theory = zeros(1,length(SNR_dB));
BER_MISO_Ala_Theory = zeros(1,length(SNR_dB));

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

figure()
semilogy(SNR_dB,BER_SISO_Theory);
hold on
semilogy(SNR_dB,BER_SIMO_MRC_Theory);
hold on
semilogy(SNR_dB,BER_MISO_MRT_Theory);
hold on
semilogy(SNR_dB,BER_MISO_Ala_Theory);
grid on;
xlabel("SNR [dB]");ylabel("Bit Error Rate [bps/Hz]");
title("Theoretical BER Result using Q-Function");
legend("SISO","SIMO","MISO MRT","MISO Alamouti")

%%
save("Theoretical_BER","BER_SISO_Theory","BER_SIMO_MRC_Theory","BER_MISO_MRT_Theory","BER_MISO_Ala_Theory");


function [BER_SISO_Theory,BER_SIMO_MRC_Theory,BER_MISO_MRT_Theory,BER_MISO_Ala_Theory]=BER_theoretical(SNR_dB_range,Nr_SIMO,Nt_MISO)
SNR_dB = 0:2:SNR_dB_range;
SNR = 10.^(SNR_dB./10);
BER_SISO_Theory = zeros(1,length(SNR_dB));
BER_SIMO_MRC_Theory = zeros(1,length(SNR_dB));
BER_MISO_MRT_Theory = zeros(1,length(SNR_dB));
BER_MISO_Ala_Theory = zeros(1,length(SNR_dB));
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