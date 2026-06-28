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
Max = 2;                            % The number fo channel realisation
a_d = 3.5;                            % path loss coefficient  for direct path
a_r = 2;                              % path loss coefficient  for reflect path
d_R = 30;                             % The distance between the transmitter and RIS (m)
d_r = 3;                              % The distance between the RIS and User (m)
dis_depe_coeff = 10^(-30/10);         % distance dependent coefficient, -30 dB to linear
kappa = 10^(3/10);                    % The Rician factor 
tolerance = 1e-8;

% Optimal
power_SOCP = zeros(Max,length(SINR_all));
power_SDP  = zeros(Max,length(SINR_all));
power_KKT  = zeros(Max,length(SINR_all));

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


     %% Effective channel
     h_k_optimal = zeros(Nt,N_UE); % with optimal phase shift design on RIS

     % Compute the effective channel for each user
     for i = 1:N_UE
        h_k_optimal(:,i) = (h_d_k(:,i)' + h_r_k(:,i)'*O_optimal*G)';
        % h_k_optimal(:,1) and h_k_random(:,1) is effective channel for UE1
     end

    % start the SINR loop
    for i_snr = 1:length(SINR_all)
       % SINR value
       SINR = SINR_all(i_snr); % set the SINR
       
        %% Precoder optimisation

        % Based ON optimised phase shift RIS design
        [power_SOCP(i_max,i_snr)]=SOCP_optimisation_RIS(h_k_optimal,SINR,sigma);     % SOCP
        [power_SDP(i_max,i_snr)]=SDP_optimisation_RIS(h_k_optimal,SINR,sigma);       % SDP
        [KKT_optimal]=KKT_optimisation_RIS(h_k_optimal,SINR,sigma,tolerance);        % KKT
        power_KKT(i_max,i_snr) = sum(KKT_optimal);

    end
end

%% Average transmit power

% optimal
Avg_p_SOCP_optimal = 10*log10(mean(power_SOCP,1))+30;
Avg_p_SDP_optimal  = 10*log10(mean(power_SDP,1))+30;
Avg_p_KKT_optimal  = 10*log10(mean(power_KKT,1))+30;

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

