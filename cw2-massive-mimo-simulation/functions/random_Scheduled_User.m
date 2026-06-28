
%?>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Input:
%
% - Scheduled_User (1xQ double) = the index of the scheduled user (1~K)
% - H_center_t(Nt*Num_User double) = The (channel) at time instant k
% - H_interfere_t(Nt*Num_User*Num_BS) = The (channel) at time instant k
% - PL_center
% - PL_interfer
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Output:
% - Scheduled_H_center_t(Nr*Nt*Q) = scheduled user channel from center BS
% - Scheduled_H_interfere_t(Nr*Nt*Q*Num_BS) = the scheduled user channel from
% interfere BS
%?>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [Scheduled_H_center_t,Scheduled_H_interfere_t]=random_Scheduled_User(Scheduled_User,H_center_t,H_interfere_t,PL_center,PL_interfer)
    [Nt,~,Num_BS] = size(H_interfere_t);
    Q = length(Scheduled_User); % determine the length
    Scheduled_H_center_t = zeros(Nt,Q);
    Scheduled_H_interfere_t = zeros(Nt,Q,Num_BS);
    
    % find the scheduled user channel
    for i_sch = 1:Q
        % center BS to Q
           Scheduled_H_center_t(:,i_sch) = H_center_t(:,Scheduled_User(i_sch))./sqrt((PL_center(Scheduled_User(i_sch))));
        
        % interfere BSs to Q
           Scheduled_H_interfere_t(:,i_sch,:) = H_interfere_t(:,Scheduled_User(i_sch),:)./sqrt((PL_interfer(:,Scheduled_User(i_sch),:)));
          
    end

end