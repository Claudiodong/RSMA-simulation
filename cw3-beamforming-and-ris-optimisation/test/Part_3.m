clc
clear
close all;
addpath("function_Part3\");
%%
orange = [0 0.4470 0.7410];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels

% Parameter
Nt = 4;                               % The number 0f transmit antenna
Nr = 1;                               % The number of receive antenna
N_UE = 4;                             % The number of users
M = 20;                               % THe number of elements in RIS
SINR_dB = -5:3:10;                    % SINR in dB
SINR_all = 10.^(SINR_dB./10);         % SINR in linear
sigma = 10.^((-70-30)/10);            % The noise power
Max = 1;                            % The number fo channel realisation
a_d = 3.5;                            % path loss coefficient  for direct path
a_r = 2;                              % path loss coefficient  for reflect path
d_R = 30;                             % The distance between the transmitter and RIS (m)
d_r = 3;                              % The distance between the RIS and User (m)
dis_depe_coeff = 10^(-30/10);         % distance dependent coefficient, -30 dB to linear
kappa = 10^(3/10);                    % The Rician factor 
tolerance = 1e-8;

% Optimal
power_SOCP_optimal = zeros(Max,length(SINR_all));
power_SDP  = zeros(Max,length(SINR_all));
power_KKT  = zeros(Max,length(SINR_all));

% Random
power_SOCP_random = zeros(Max,length(SINR_all));
power_SDP_random  = zeros(Max,length(SINR_all));
power_KKT_random  = zeros(Max,length(SINR_all));

% MU MISO
SOCP_MU_MISO = zeros(Max,length(SINR_all));
SDP_MU_MISO = zeros(Max,length(SINR_all));
KKT_MU_MISO = zeros(Max,length(SINR_all));

% Start the realisation loop
i_max = 0;
while(i_max<Max)
    i_max = i_max +1 ; % counter
    sprintf("Running the %1.1f Realisation",i_max)

     % Direct and Reflect path
     [h_d_k,h_r_k]=Channel_d_r(Nt,N_UE,a_d,d_r,d_R,dis_depe_coeff,a_r,kappa,M);

     % Channel between the transmitter to RIS
     G=G_channel(Nt,N_UE,M,dis_depe_coeff,d_R,a_r,kappa);

     %% phase shift optimisation
     [O_optimal]=RIS_phase_optimisation(h_d_k,h_r_k,G);

     % Random phase shift generation    
     theta = (2*pi)*rand(M,1);  % random generate the theta value from [0,2*pi]
     O_random = diag(exp(1i*theta));% convert to complex form

     %% Effective channel
     h_k_optimal = zeros(Nt,N_UE); % with optimal phase shift design on RIS
     h_k_random = zeros(Nt,N_UE);  % Random phase shift design on RIS

     % Compute the effective channel for each user
     for i = 1:N_UE
        h_k_optimal(:,i) = (h_d_k(:,i)' + h_r_k(:,i)'*O_optimal*G)';
        % h_k_optimal(:,1) and h_k_random(:,1) is effective channel for UE1
        h_k_random(:,i) =  (h_d_k(:,i)' + h_r_k(:,i)'*O_random*G)';          
     end

    % start the SINR loop
    for i_snr = 1:length(SINR_all)
       % SINR value
       SINR = SINR_all(i_snr); % set the SINR
       
        %% Precoder optimisation

        % Based ON optimised phase shift RIS design
        [power_SOCP_optimal(i_max,i_snr)]=SOCP_optimisation_RIS(h_k_optimal,SINR,sigma);     % SOCP
        [power_SDP(i_max,i_snr)]=SDP_optimisation_RIS(h_k_optimal,SINR,sigma);       % SDP
        [KKT_optimal]=KKT_optimisation_RIS(h_k_optimal,SINR,sigma,tolerance);        % KKT
        power_KKT(i_max,i_snr) = sum(KKT_optimal);

        %% Random phase shift
        [power_SOCP_random(i_max,i_snr)]=SOCP_optimisation_RIS(h_k_random,SINR,sigma);    
        [power_SDP_random(i_max,i_snr)]=SDP_optimisation_RIS(h_k_random,SINR,sigma);       % SDP
        [KKT_randoml]=KKT_optimisation_RIS(h_k_random,SINR,sigma,tolerance);        % KKT
        power_KKT_random(i_max,i_snr) = sum(KKT_randoml);

        %% MU MISO (checked)
        [SOCP_MU_MISO(i_max,i_snr)]=SOCP_optimisation_MU_MISO(h_d_k,sigma,SINR);% SOCP
        [SDP_MU_MISO(i_max,i_snr)]=SDP_optimisation_RIS(h_d_k,SINR,sigma);      % SDP    
        [KKT_MU]=KKT_optimisation_RIS(h_d_k,SINR,sigma,tolerance);              % KKT
        KKT_MU_MISO(i_max,i_snr) = sum(KKT_MU);    
    end
end

%% Average transmit power

% optimal
Avg_p_SOCP_optimal = 10*log10(mean(power_SOCP_optimal,1))+30;
Avg_p_SDP_optimal  = 10*log10(mean(power_SDP,1))+30;
Avg_p_KKT_optimal  = 10*log10(mean(power_KKT,1))+30;

% Random
Avg_p_SOCP_random = 10*log10(mean(power_SOCP_random,1))+30;
Avg_p_SDP_random  = 10*log10(mean(power_SDP_random,1))+30;
Avg_p_KKT_random  = 10*log10(mean(power_KKT_random,1))+30;

% MU MISO case without the RIS (checked correct)
Avg_p_SOCP_MU_MISO = 10*log10(mean(SOCP_MU_MISO,1))+30;
Avg_p_SDP_MU_MISO  = 10*log10(mean(SDP_MU_MISO,1))+30;
Avg_p_KKT_MU_MISO  = 10*log10(mean(KKT_MU_MISO,1))+30;

%% plot the figures
figure()
% dBm plot
plot(SINR_dB,Avg_p_SOCP_optimal,"k-square",'LineWidth',1.5,'MarkerSize',10)
hold on;
plot(SINR_dB,Avg_p_SDP_optimal,"r--O",'LineWidth',1.5,'MarkerSize',10)
hold on;
plot(SINR_dB,Avg_p_KKT_optimal,"-.*","Color",orange,'LineWidth',1.5,'MarkerSize',10)
hold on;
grid on;xlabel("SINR Requirment [dB]");ylabel("Transmit Power [dBm]")
legend("SOCP","SDP","KKT",'Location','northwest')

%% Comparsion
figure()
% watt plot
plot(SINR_dB,Avg_p_SOCP_optimal,"r-square",'LineWidth',1.5,'MarkerSize',13)
hold on;
plot(SINR_dB,Avg_p_SDP_optimal,"g-O",'LineWidth',1.5,'MarkerSize',8)
hold on;
plot(SINR_dB,Avg_p_KKT_optimal,"b-x",'LineWidth',1.5,'MarkerSize',8)
hold on;
plot(SINR_dB,Avg_p_SOCP_random,"--square",'color',orange,'LineWidth',1.5,'MarkerSize',13)
hold on;
plot(SINR_dB,Avg_p_SDP_random,"m--O",'LineWidth',1.5,'MarkerSize',8)
hold on;
plot(SINR_dB,Avg_p_KKT_random,"k--x",'LineWidth',1.5,'MarkerSize',8)
hold on;
plot(SINR_dB,Avg_p_SOCP_MU_MISO,"-.square",'LineWidth',1.5,'MarkerSize',13)
hold on;
plot(SINR_dB,Avg_p_SOCP_MU_MISO,"-.O",'LineWidth',1.5,'MarkerSize',8)
hold on;
plot(SINR_dB,Avg_p_SOCP_MU_MISO,"-.x",'LineWidth',1.5,'MarkerSize',8)
grid on;xlabel("SINR Requirment [dB]");ylabel("Transmit Power [dBm]");
legend("SOCP Optimal","SDP Optimal",'KKT Optimal','SOCP Random','SDP Random','KKT Random','SOCP MU-MISO','SDP MU-MISO','KKT MU-MISO')


