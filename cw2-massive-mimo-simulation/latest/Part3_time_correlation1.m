clc
clear
close all;
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
shadowing_correlation = 0;
Nt = 4;                                 % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr1 = [0,0.5,0.98];                % time correlation for channels (time correlation)
time_k = 1e3;                             % time instant (spatial correlation)
drop = 1e3;                             % Number of drop
Avg_C_time = zeros(drop,Num_User);
Avg_C_corr = zeros(length(t_corr1),drop*Num_User);

for i_c = 1:length(t_corr1)
    t_corr = t_corr1(i_c);
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
             Q = 4;                                  % the number of scheduled user at this time instant
             
             % The index of the scheduled User
             Scheduled_User = randperm(Num_User,Q);
    
             Scheduled_H = zeros(Nr,Nt,Q,Num_BS+1);
             % scheduled user channel
             for i_q = 1:Q
                 for i_bs = 1:Num_BS+1
                     Scheduled_H(:,:,i_q,i_bs) = H(:,:,Scheduled_User(i_q),i_bs)/sqrt(PL(i_bs,Scheduled_User(i_q)));
                 end
    
             end

    
             % Compute the ZFBF for the center BS
             center_H = squeeze(Scheduled_H(:,:,:,1)).'; % need to transpose into Q*Nt
             F = center_H'*inv(center_H*center_H');
             w = zeros(Nt,Q);
             for i_Q = 1:Q        
                 % precoder design for the center BS
                 w(:,i_Q) = F(:,i_Q)/norm(F(:,i_Q));
             end
    
             % Uniform Power Allocation 
             Ne = min(min(Nt,Num_User),Q);
             S = Tx_power/Ne;
   
             % Random precoder design for the interfere BSs
             H_interfere = pagetranspose(squeeze(Scheduled_H(:,:,:,2:end))); % transpose into Q*Nt*Num_BS
             w_interfere = (randn(Nt,Q,Num_BS) + 1i*randn(Nt,Q,Num_BS))/sqrt(2) ;
             A = zeros(Q,Num_BS);
             for i_bs = 1:Num_BS
                 A(:,i_bs) = diag(S.*abs(H_interfere(:,:,i_bs)*w_interfere(:,:,i_bs)).^2);
             end
    
    
             % Compute the SINR for the User
             UE_all_power = S*abs(center_H*w).^2;
             UE_power = diag(UE_all_power);
    
             for i_e = 1:Q
                 % inter-user interference (Co-scheduled), should be close to 0
                 % by ZFBF
                 I_c = sum(UE_all_power(i_e,:)) - UE_power(i_e);
    
                 % interference from all 6 interfere BS to UE
                 I_interfere = sum(A(i_e,:));
    
                 % compute the SINR
                 SINR = UE_power(i_e)/(N0+I_c+I_interfere);
    
                 % Compute the Rate
                 C = log2(1+SINR);
    
                 % Mapping back to the corresponding user
                 Sum_C(i_t,Scheduled_User(i_e)) = C;
             end
             
             
             % Update the Small scale fading
             N = (randn(Nr,Nt,Num_User,Num_BS+1) + 1i* randn(Nr,Nt,Num_User,Num_BS+1))/sqrt(2);
             H = t_corr.*H + sqrt(1-t_corr^2).*N;
             
         end
         % Average User Rate per drop after time_k time
         Avg_C_time(i_drop,:) = sum(Sum_C)./time_k;
    end
    Avg_sort = sort(reshape(Avg_C_time,1,[]),"ascend");
    Avg_C_corr(i_c,:) = Avg_sort;
end
%% plot

CDF = (1:length(Avg_sort))./length(Avg_sort) *100;

figure()
plot(Avg_C_corr,CDF,"--","LineWidth",1.5)
grid on;xlabel("Average User Rate [bps/Hz]");ylabel("Cumulative Distribution Function [\%]");
legend("$\epsilon=0$", "$\epsilon=0.85$", "$\epsilon=0.98$")
title("")

