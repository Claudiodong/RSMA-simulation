
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Input:
% - Nr(1x1 double)= The number of receive antenna on UEs
%
% - Nt(1x1 double) = The number transmit antenna on each BS
%
% - Num_User(1x1 double) = the number of user that require in the system
% - Num_BS(1x1 double) = the number of interfer BS in the system
% - time_k(1x1 double) = the time instant, how long the drop is
% - t_corr(1x1 double) = the time correlation of channel direction
% information from the previou to the present H^hat(t) =
% u*H^hat(t-1)+sqrt(1-u^2)*N(t)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%
% - H_center (Nr*Nt*Num_User*time_k) = the channel for the user from
% the center BS
%
% - H_interfer（Nr*Nt*Num_User*Num_BS*time_k）= the channel for each user from
% each BS for the all time_k 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [H_center,H_interfer] = small_scale_channel(Nr,Nt,Num_User,Num_BS,time_k,t_corr)

    %time_k = 4; % the number of time step require
    % NrxNt x Number of User x Number of BS x time => 
    % (:,:,2,1,1) = the channel matrix for the first interfer BS to the second
    % user at time instant k = 1.
    N_interfer = zeros(Nr,Nt,Num_User,Num_BS,time_k);
    N_center = zeros(Nr,Nt,Num_User,time_k);
    
    H_interfer = zeros(Nr,Nt,Num_User,Num_BS,time_k);
    H_center = zeros(Nr,Nt,Num_User,time_k);

    % At time instant k = 1, since matlab index count from 1. based on
    % randn
    N_interfer(:,:,:,:,1) = 1/sqrt(2) * (randn(Nr,Nt,Num_User,Num_BS)+1i*randn(Nr,Nt,Num_User,Num_BS));
    N_center(:,:,:,1) =     1/sqrt(2) * (randn(Nr,Nt,Num_User)+1i*randn(Nr,Nt,Num_User));

    H_interfer(:,:,:,:,1) = 1/sqrt(2) * (randn(Nr,Nt,Num_User,Num_BS)+1i*randn(Nr,Nt,Num_User,Num_BS));
    H_center(:,:,:,1) =     1/sqrt(2) * (randn(Nr,Nt,Num_User)+1i*randn(Nr,Nt,Num_User));
        
    for i = 2:time_k
        % The new N
        N_interfer(:,:,:,:,i) = 1/sqrt(2) * (randn(Nr,Nt,Num_User,Num_BS)+1i*randn(Nr,Nt,Num_User,Num_BS));
        N_center(:,:,:,i) =     1/sqrt(2) * (randn(Nr,Nt,Num_User)+1i*randn(Nr,Nt,Num_User));
    
        % The new channel infromation
        H_interfer(:,:,:,:,i) = t_corr*H_interfer(:,:,:,:,i-1) + sqrt(1-(t_corr^2))*N_interfer(:,:,:,:,i);
        H_center(:,:,:,i) =     t_corr*H_center(:,:,:,i-1) + sqrt(1-(t_corr^2))*N_center(:,:,:,i);
    end
end