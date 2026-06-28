clc
clear
close all;
orange = [0.8500 0.3250 0.0980];
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((N0_dBm-30)/10);               % noise variance in linear
Num_User = 10;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
shadowing_correlation = 0;
Nt = 4;                                % Number of Tx antenna
Nr = 2;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels (time correlation)
time_k = 3;                             % time instant (spatial correlation)
drop = 1e3;                             % Number of drop
Q = 3;                                  % Number of Scheduled User

%
avg_C = zeros(drop,Q);

    for i_drop = 1:drop  % Number of drop loop, where the UEs are at fixed location
    
        %Random deploy User
    
        % center BS to users
        [User_location,Center_BS_UEs] = Random_user_deploy(d(2),d(1),Num_User,center);
            
        % interfer BS location
        BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
        BS_location = zeros(2,Num_BS);
        for i = 1:Num_BS
             BS_location(:,i) = r*[cos(BS_theta(i));sin(BS_theta(i))];
        end
    
        %% Large Scale Fading do not change since the UEs are fixed within a Drop time
            
        % center BS to Users
        % Path loss large scale fading in km and dB, Center BS to UE
        PL_dB = 128.1+37.6*log10(Center_BS_UEs./1000)+shadowing_dB*randn(Nr,Num_User);
        % Path loss in linear, Center BS to UEs
        PL_center_BS_user = 10.^((PL_dB)./10);
            
        % other BS to users Path loss
        interfer_BS_to_UE = zeros(1,Num_User,Num_BS);
        for i = 1:Num_BS
            interfer_BS_to_UEs = BS_location(:,i) - User_location; % Every BS disatnce to all users,[BS1-BS6]
            interfer_BS_to_UE(:,:,i) = sqrt(sum(interfer_BS_to_UEs.^2)); % distance
        end
            
        % Path loss in dB
        PL_dB_other_BS_uEs = 128.1+37.6*log10(interfer_BS_to_UE./1000)+shadowing_dB*randn(Nr,Num_User,Num_BS);
        % Path loss in linear, other BS to UEs
        PL_interfer_BS_UE = 10.^((PL_dB_other_BS_uEs)./10);    
        %PL_interfer_BS_UE = squeeze_size(PL_interfer_BS_UE); % squeeze the size if there is diemension of 1
        % PL_interfer_BS_UE = Num_UEs * Num_BS
        % therefore, (:,1) = The Interfere BS 1 to all UEs
        % (1,:) = The All path loss from all interfere BS to UE1.
            
        %% small scale fading channel
        [H_center,H_interfer] = small_scale_channel(Nr,Nt,Num_User,Num_BS,time_k,t_corr);
        % H_interfer is a 5-D matrix, where H_interfer = NrxNr x User x Num_BS x Time instant
        % H_center => (Nr*Nt*Num_User*time_k)
    
        C = zeros(1,Q);
        % spatial correlation
        for i_time = 1:time_k % the number of time instant that the channel will vary when the UEs are fixed
             H_center_t = H_center(:,:,:,i_time);
             H_interfere_t = H_interfer(:,:,:,:,i_time);
             % scheduled User channel
             Scheduled_User = randperm(Num_User,Q);

             % fine the corresponding user channel from center and
             % interfere BSs
             %  have size of Nt*Nr*Q*(Num_BS)
             [Scheduled_H_center_t,Scheduled_H_interfere_t]=random_Scheduled_User(Scheduled_User,H_center_t,H_interfere_t,PL_center_BS_user,PL_interfer_BS_UE);
              
             % Transpose into Nr*Nt*Q*Num_BS form 
             Scheduled_H_center_t1 = pagetranspose(Scheduled_H_center_t);
             Scheduled_H_interfere_t1 = pagetranspose(Scheduled_H_interfere_t);

             % block diagonal 
             [U,S,V] = pagesvd(Scheduled_H_center_t1); % checked, works
             % find the null subspace
             Anull_subapce = V(:,Nr+1:Nt,:);

             for i = 1:Q
                 A = Scheduled_H_center_t1(:,:,i)*Anull_subapce(:,:,i);
                 [~,~,V1] = pagesvd(A);
                 V_domaint = V1(:,1);
                 % Nt*Q
                 % precoder
                 W = Anull_subapce(:,:,i)*V_domaint;
             end
                        
        end
        
    end

function [Scheduled_H_center,Scheduled_H_interfere]=random_Scheduled_User(Scheduled_User,H_center_t,H_interfer_t,PL_center_BS_user,PL_interfer_BS_UE)
             % find the parameters
             [Nr,Nt,~,Num_BS] = size(H_interfer_t);

             % define the size
             Q = length(Scheduled_User);

             % obtain the scheduled user interference from iinterfere BS
             Scheduled_H_interfere = zeros(Nt,Nr,Q,Num_BS);
             Scheduled_H_center =  zeros(Nt,Nr,Q);

             for i_UE = 1:Q
                 Scheduled_U = Scheduled_User(i_UE);

                 A = pagetranspose(H_interfer_t(:,:,Scheduled_U));
                 B = PL_interfer_BS_UE(:,Scheduled_U,:);

                 % Interfere BS channels
                 Scheduled_H_interfere(:,:,i_UE,:) = A./sqrt(B);

                 % Center BS for scheduled user
                 Scheduled_H_center(:,:,i_UE) = (H_center_t(:,:,Scheduled_U).')./sqrt(PL_center_BS_user(:,Scheduled_U));
             end
end

function [w_center,w_interfere]=ZFBF_MIMO(Scheduled_H_center_t,Scheduled_H_interfere_t)

            % Define the variables
            [Nr,Nt,Q,Num_BS] = size(Scheduled_H_interfere_t);
    
%             F_center = zeros(Nt,Nr,Q);
%             F_interfere = zeros(Nt,Nr,Q,Num_BS);

            % Zero-forcing beamforming, w for the center BS to scheduled
            % user
            w_center = zeros(Nt,Nr,Q);
            % find the interfere F
             for i_UE = 1:Q
                 F_center= pinv(Scheduled_H_center_t(:,:,i_UE));
%                  for i_sc = 1:Num_BS
%                      % ZFBF for each user from the interfere BS
%                      F_interfere(:,:,i_UE,i_sc) = pinv(Scheduled_H_interfere_t(:,:,i_UE,i_sc));
%                  end
                 w_center(:,:,i_UE) = F_center./norm(F_center);
                 
             end
                 
             % Random precoding for the interfere BS to users
             w_interfere = (randn(Nt,Nr,Q,Num_BS) + 1i*randn(Nt,Nr,Q,Num_BS))/sqrt(2);

end