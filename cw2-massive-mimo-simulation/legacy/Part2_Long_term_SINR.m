clc
clear
close all;
orange = [0.8500 0.3250 0.0980];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((-174-30)/10);                 % noise variance in linear
Num_User = 1e4;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
shadowing_correlation = 0;
Nt = 4;                                 % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels
time_k = 4;                             % time instant
MAX = [1,10,100];                       % Number of Monte Carlo Runs
avg_dB = zeros(length(MAX),Num_User);
cdf = zeros(length(MAX),Num_User);

for i_max = 1:length(MAX)
    i_runs = 0;
    sum_SINR = 0;
    % Monte Carlo
    while(i_runs<MAX(i_max))
        i_runs = i_runs +1;
    
        %% Random deploy User
        % center BS to users
        [User_location,Center_BS_UEs] = Random_user_deploy(d(2),d(1),Num_User,center);
        
        % interfer BS location
        BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
        BS_location = zeros(2,1,Num_BS);
        for i = 1:Num_BS
            BS_location(:,:,i) = r*[cos(BS_theta(i));sin(BS_theta(i))];  
        end
        
        %% Large Scale Fading
        
        % center BS to Users
        % Path loss large scale fading in km and dB, Center BS to UE
        PL_dB = 128.1+37.6*log10(Center_BS_UEs/1000)+shadowing_dB*randn(1,Num_User);
        % Path loss in linear, Center BS to UEs
        PL_center_BS_user = 10.^((PL_dB)./10);
        
        % other BS to users Path loss 
        interfer_BS_to_UEs = BS_location - User_location; % Every BS disatnce to all users,[BS1-BS6]
        interfer_BS_to_UEs = sqrt(sum(interfer_BS_to_UEs.^2)); % distance
        
        % Path loss in dB
        PL_dB_other_BS_uEs = 128.1+37.6*log10(interfer_BS_to_UEs./1000)+shadowing_dB*randn(1,Num_User,Num_BS);
        % Path loss in linear, other BS to UEs
        PL_interfer_BS_UE = 10.^((PL_dB_other_BS_uEs)./10);
        
        % PL_interfer_BS_UE => Nt x Number of User x Number of BS
        % e.g (:,:,1) => the all path loss for the first BS to all 10 users
        
        %% Long-term SINR
        
        % Uniform power allocation
        E_UE_H = Tx_power;   % Transmit power per UE 
        %E_UE_H = Tx_power/Nt;   % Transmit power per UE 
        
        % for Q2 => ONLY consider the Large sacle fading
        H_largesacle_center_BS = 1./PL_center_BS_user;
        H_largescale_interfer_BS = 1./PL_interfer_BS_UE;
        
        SINR_Center_BS_UEs = (H_largesacle_center_BS.*E_UE_H)./(N0 + sum(H_largescale_interfer_BS.*E_UE_H,3));
        % sum(H_largescale_interfer_BS.*E_UE_H,3) compute a result of 1xNum_UE,
        % since it need to find all the interference to this User from all the BSs
        % eg. for UE1, it have interference from all interfer BS and noise, which
        % experience different path loss
        
        %SINR_Center_BS_UE1 = (H_largesacle_center_BS(1).*E_UE_H)./(N0 + sum(H_largescale_interfer_BS(:,1,:).*E_UE_H));
        % using this to check the values
    
        sum_SINR = SINR_Center_BS_UEs + sum_SINR;
    end
    % average SINR and sorted from low to high
    avg_dB(i_max,:) = sort(10*log10(sum_SINR./MAX(i_max)));
    % compute the CDF
    cdf(i_max,:) = (1:length(avg_dB))./(length(avg_dB))*100;
end
%%
figure()
plot(avg_dB(1,:),cdf(1,:),"-","LineWidth",3);grid on;
hold on;plot(avg_dB(2,:),cdf(2,:),"-","LineWidth",3);
hold on;plot(avg_dB(3,:),cdf(3,:),"-","LineWidth",3);
xlabel("Long-Term SINR at Users [dB]");ylabel("Cumulative Distribution Function [\%]");
legend("Monte Carlo = 1","Monte Carlo = 10","Monte Carlo = 100","Location","northwest");
title("The Long-Term SINR CDF with differet Monte Carlo Number")


