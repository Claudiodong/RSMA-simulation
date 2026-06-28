clc
clear
close all;
addpath("New Function\");
orange = [0.8500 0.3250 0.0980];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((N0_dBm-30)/10);               % noise variance in linear
Num_User = 10;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
Nt1 = [4,16,64];                        % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels (time correlation)
time_k = 10;                           % time instant (spatial correlation)
drop = 1e3;                             % Number of drop
Q = 4;                                  % the number of scheduled user at this time instant
Avg_C_Nt = zeros(length(Nt1),drop*Num_User);
Avg_C_time = zeros(drop,Num_User);

for i_nt= 1:length(Nt1)
    Nt = Nt1(i_nt);
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
         PL = 10.^((PL_dB)./10);
         
         % initilise the small scale channel 
         %[~,H] = small_scale_channel(Nr,Nt,Num_User,Num_BS+1,time_k,t_corr);
         % Time instant loop where the small scale fading change
         H = (randn(Nr,Nt,Num_User,Num_BS+1) + 1i* randn(Nr,Nt,Num_User,Num_BS+1))/sqrt(2);
         
         Sum_C = zeros(time_k,Num_User);
         for i_t = 1:time_k
             % Random User Schedule
             
             % The index of the scheduled User
             Scheduled_User = randperm(Num_User,Q);
    
             % Scheduled User Channels for all BS(including the center BS)
             [Scheduled_H]=Scheduled_User_H(Scheduled_User,H,PL);

             % transmit power
             S = Tx_power/Q;
    
             % Compute the ZFBF
             center_H = squeeze(Scheduled_H(:,:,:,1)).'; % need to transpose into Q*Nt
             H_interfere = pagetranspose(squeeze(Scheduled_H(:,:,:,2:end))); % transpose into Q*Nt*Num_BS
    
             % The ZFBF w matrix
             [w_center,w_interfere]=ZFBF_scheduled_UEs(center_H,H_interfere);
    
             % capacity 
             C = Capacity(center_H,H_interfere,w_center,w_interfere,S,N0);
             % Mapping back to corresponding UE capacity location
             for i= 1 : Q
                 Sum_C(i_t,Scheduled_User(i)) = C(i);
             end
             
             % Update the Small scale fading
             [H]=Update_Small_Scale_Fading(H,t_corr);        
         end
         % Average User Rate per drop after time_k time
         Avg_C_time(i_drop,:) = sum(Sum_C,1)./time_k;
    end
    Avg_C_Nt(i_nt,:) = sort(reshape(Avg_C_time,1,drop*Num_User),"ascend");
end

%% plot
CDF = (1:length(Avg_C_Nt))./length(Avg_C_Nt) *100;

figure()
for i = 1:length(Nt1)
plot(Avg_C_Nt(i,:),CDF)
hold on;
end
grid on;
xlabel("Average User Rate [bps/Hz]");ylabel("Cumulative Distribution Function [\%]")
legend

save("Nr1Q4","Avg_C_Nt")