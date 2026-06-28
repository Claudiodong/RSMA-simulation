clc
clear
close all;
addpath("function\")
%% parameter
Tx_power = 10^((46-30)/10);             % Transmit power in linear
N0 = 10^((-174-30)/10);                 % noise variance in linear
Num_User = 10;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
Nt1 = [4,16,64];                        % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels (time correlation)
time_k = 10;                            % time instant (spatial correlation)
drop = 1e3;                             % Number of drop
Avg_C = zeros(drop,Num_User);
Avg_C_all = zeros(length(Nt1),drop*Num_User);

for i_Nt = 1:length(Nt1)
    for i_drop = 1:drop                 % Number of drop loop, where the UEs are at fixed location
        %% Random deploy User
        Nt = Nt1(i_Nt);
        % center BS to users
        [User_location,Center_BS_UEs] = Random_user_deploy(d(2),d(1),Num_User,center);
            
        % interfer BS location
        BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
        BS_location = zeros(2,Nr,Num_BS);
        for i = 1:Num_BS
             BS_location(:,:,i) = r*[cos(BS_theta(i));sin(BS_theta(i))];
        end
    
        %% Large Scale Fading do not change since the UEs are fixed within a Drop time
            
        % center BS to Users
        % Path loss large scale fading in km and dB, Center BS to UE
        PL_dB = 128.1+37.6*log10(Center_BS_UEs./1000)+shadowing_dB*randn(Nt,Num_User);
        % Path loss in linear, Center BS to UEs
        PL_center_BS_user = 10.^((PL_dB)./10);
            
        % other BS to users Path loss 
        interfer_BS_to_UEs = BS_location - User_location; % Every BS disatnce to all users,[BS1-BS6]
        interfer_BS_to_UEs = sqrt(sum(interfer_BS_to_UEs.^2)); % distance
            
        % Path loss in dB
        PL_dB_other_BS_uEs = 128.1+37.6*log10(interfer_BS_to_UEs./1000)+shadowing_dB*randn(Nt,Num_User,Num_BS);
        % Path loss in linear, other BS to UEs
        PL_interfer_BS_UE = 10.^((PL_dB_other_BS_uEs)./10);    

            
        %% Small scale fading channel
        [H_center,H_interfer] = small_scale_channel(Nr,Nt,Num_User,Num_BS,time_k,t_corr);
        % H_interfer is a 5-D matrix, where H_interfer = NrxNr x User x Num_BS x Time instant
        % H_center => (Nr*Nt*Num_User*time_k)
        C_total = zeros(1,Num_User);
        C_time = zeros(time_k,Num_User);
        % spatial correlation
        for i_time = 1:time_k % the number of time instant that the channel will vary when the UEs are fixed
             H_center_time = squeeze_size(H_center(:,:,:,i_time)); % have size Nt x Num_UEs
             H_interfer_time = squeeze_size(H_interfer(:,:,:,:,i_time));% have size Nt x Num_UEs x Num_BS

             Q = 2;                  % Number of Scheduled User
              
             %% Random scheduled User channel
             Scheduled_User = randperm(Num_User,Q);
         
             %% find the scheduled user channels
             [Scheduled_H_center_t,Scheduled_H_interfere_t]=random_Scheduled_User(Scheduled_User,H_center_time,H_interfer_time,PL_center_BS_user,PL_interfer_BS_UE);
    
             % Transpose into Q*Nt form 
             Scheduled_H_center_t = Scheduled_H_center_t.';
             Scheduled_H_interfere_t = pagetranspose(Scheduled_H_interfere_t);
          
             % Perform transmit schemes at all BS (Zero-Forcing Beamforming ZFBF)
             % the precoding scheme for the interfere BS is random precoding
             [w_center,w_interfere]=ZFBF_sechulded_user(Scheduled_H_center_t,Scheduled_H_interfere_t);
                        
             % Compute the SINR for the scheduled users
    
             % Uniform power allocation for the scheduled users, since 1x4 set
             % up => SIMO system, not inter-stream interference between the
             % stream
             Ne = min(min(Nt,Num_User),Q);
             Es_per_channel = (Tx_power/Q); % transmit power per UE per channel

             [SINR]=instantaneous_SINR(Scheduled_H_center_t,Scheduled_H_interfere_t,w_center,w_interfere,Es_per_channel,N0);
    
             % Rate and accumulate it
             % scheduled user rate
             C = log2(1+SINR);

             % Map into each user
             for i_q = 1:Q
               C_total(Scheduled_User(i_q)) = C_total(Scheduled_User(i_q))+C(i_q);
               C_time(i_time,Scheduled_User(i_q)) = C(i_q);
             end
             
        end
        Avg_C(i_drop,:) = C_total./time_k;
    end
   Avg_C_all(i_Nt,:)=sort(reshape(Avg_C,1,[]));
end
%%
cdf = (1:length(Avg_C_all))./length(Avg_C_all)*100;

figure()
plot(Avg_C_all,cdf)
grid on;
ylabel("CDF [%]")
xlabel("Average User Rate [bps/Hz]")
legend("Nt=4","Nt=16","Nt=64")

