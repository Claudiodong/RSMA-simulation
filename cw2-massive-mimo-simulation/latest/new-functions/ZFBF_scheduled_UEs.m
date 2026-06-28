
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - center_H (Q*Nt)= The center BS channel to scheduled users 
% - H_interfere (Q*Nt*Num_BS)= The center BS channel to scheduled users 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - w_center(Nt*Q) = The w ZFBF precoder design for the center BS to all
% scheduled Users
%
% - w_interfere(Nt*Q*Num_BS) = The w ZFBF precoder design for the interfere BS to all
% scheduled Users
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [w_center,w_interfere]=ZFBF_scheduled_UEs(center_H,H_interfere)

         % Find the Parameters
         [Q,Nt,Num_BS]=size(H_interfere); % Here the Num_BS is 6, only the interfere BSs

         % Compute the ZFBF for the center BS
         F = center_H'*inv(center_H*center_H');
         w_center = zeros(Nt,Q);
         for i_Q = 1:Q        
             % precoder design for the center BS
             w_center(:,i_Q) = F(:,i_Q)/norm(F(:,i_Q));
         end
         
         % Random precoder design for the interfere BSs
         w_interfere = (randn(Nt,Q,Num_BS) + 1i*randn(Nt,Q,Num_BS))/sqrt(2) ;
         % normalise the random precoder
         w_interfere = w_interfere./pagenorm(w_interfere);

end