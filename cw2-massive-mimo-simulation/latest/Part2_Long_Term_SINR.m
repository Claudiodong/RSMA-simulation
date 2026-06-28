clc
clear
close all;
addpath("New Function\");

% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((N0_dBm-30)/10);               % noise variance in linear
Num_User = 1e5;                         % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
shadowing_correlation = 0;
Nt = 4;                                 % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels (time correlation)
time_k = 10;                            % time instant (spatial correlation)
drop = 1;                               % Number of drop
Q = 4;                                  % the number of scheduled user at this time instant


% drops loop, where the UEs locations unchange,


    for i_drop = 1:drop    
         % center BS to users
         [User_location,center_BS] = Random_user_deploy(d(2),d(1),Num_User,center);
                
         % interfer BS location
         BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
         BS_location = zeros(2,Num_BS+1);
         Distance = zeros(Num_BS+1,Num_User);
    
         for i = 2:Num_BS+1
              % interfere BSs location 
              BS_location(:,i) = r*[cos(BS_theta(i-1));sin(BS_theta(i-1))];
    
              % distance 
              Distance(i,:) = sqrt(sum((BS_location(:,i) - User_location).^2));
         end
         Distance(1,:) = center_BS;
    
         % Compute the path loss in dB
         PL_dB = 128.1+37.6*log10(Distance./1000)+shadowing_dB*randn(Num_BS+1,Num_User);
         % Path loss in linear, Center BS to UEs
         PL = 10.^(PL_dB./10);
    
         E_UE_H = Tx_power/Num_User;


         % interference power
         H_inter = sum((E_UE_H./PL(2:end,:)),1);
    
         % Long-Term SINR of center BS to UEs        
         L_SINR = (E_UE_H./PL(1,:))./(N0 + H_inter);

         L_SINR_dB = 10.*log10(L_SINR);
             
    end
    Avg_SINR = sort(L_SINR_dB,"ascend");
    CDF = (1:length(Avg_SINR))./length(Avg_SINR)*100;
%%
LW = 3;
figure()
plot(Avg_SINR,CDF,"LineWidth",LW);grid on;
hold on;
xlabel("Long-Term SINR of Users [dB]");ylabel("Cumulative Distribution Function [\%]");
%title("The CDF of User Long-Term SINR")



