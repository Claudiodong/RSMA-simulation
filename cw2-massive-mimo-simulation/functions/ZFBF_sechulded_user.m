%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Performing Zero-forcing Beamforming (ZFBF) for the sechulded user(Q)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Input:
% - Scheduled_H_center_t(Q*Nt) = scheduled user channel at time
% instant k from the center BS
%
% - Scheduled_H_interfere_t(Q*Nt*Num_BS) = scheduled user channel at time
% instant k from all interfere BSs
%
% where the Q is the number of scheduled user for the transmission at time
% instant k
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
% - w_center(Nt*Q) = ZFBF precoder at center BS for all Q users
% - w_interfere(Nt*Q*Num_BS) = ZFBF precoder design at interfere BSs for
% all Q users 
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [w_center,w_interfere]=ZFBF_sechulded_user(Scheduled_H_center_t,Scheduled_H_interfere_t)

        % Define the variables
        [Q,Nt,Num_BS] = size(Scheduled_H_interfere_t);

        F_center = pinv(Scheduled_H_center_t);
        F_interfere = zeros(Nt,Q,Num_BS);
        % find the interfere F
         for i_sc = 1:Num_BS
             % ZFBF for each user from the interfere BS
             F_interfere(:,:,i_sc) = pinv(Scheduled_H_interfere_t(:,:,i_sc));
             % (:,:,1) = interfere BS 1 to all four users ZFBF
         end

         w_center = zeros(Nt,Q);
         
         % Random precoding for the interfere BS to users
         w_interfere = (randn(Nt,Q,Num_BS) + 1i*randn(Nt,Q,Num_BS))/sqrt(2) ;
         for i=1:Q
            w_center(:,i) = F_center(:,i)./norm(F_center(:,i));
%             for j = 1:Num_BS
%                 w_interfere(:,i,j) = F_interfere(:,i,j)./norm(F_interfere(:,i,j));
%             end
         end
end