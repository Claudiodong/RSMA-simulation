
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Scheduled_User(1xQ Double) = The index of the scheduled user
%
% - H (Nr*Nt*Num_UEs*Num_BS Complex) = The samll scale fading for all BS,
%   (The first (:,:,:,1) is center BS)
%
% - PL (Num_BS*Num_UEs Double) = The (Large scale fading)
% Path Loss for all BS including the center BS for all Users.
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - Scheduled_H (Nr*Nt*Q*Num_BS Complex) = The Scheduled User Channel
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [Scheduled_H]=Scheduled_User_H(Scheduled_User,H,PL)

         % Here the Num_BS include the center BS,(total 7 BS)
         [Nr,Nt,~,Num_BS]=size(H);
         % Number of Scheduled User
         Q = length(Scheduled_User);

         Scheduled_H = zeros(Nr,Nt,Q,Num_BS);
         % Scheduled user channel
         for i_q = 1:Q
             for i_bs = 1:Num_BS
                 Scheduled_H(:,:,i_q,i_bs) = H(:,:,Scheduled_User(i_q),i_bs)/sqrt(PL(i_bs,Scheduled_User(i_q)));
             end

         end

end