
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - H (Nr*Nt*Num_UEs*Num_BS Complex) = The samll scale fading for all BS,
%   (The first (:,:,:,1) is center BS)
%
% - t_corr(1x1 double) =  The time correlation between the new and past
% small scale fading in time.
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - H (Nr*Nt*Num_UEs*Num_BS Complex) = The new time samll scale fading for all BS,
%   (The first (:,:,:,1) is center BS)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [H]=Update_Small_Scale_Fading(H,t_corr)

         % The Num_BS including the center BS (7 BS)
         [Nr,Nt,Num_User,Num_BS]=size(H);
         % Update the Small scale fading
         N = (randn(Nr,Nt,Num_User,Num_BS) + 1i* randn(Nr,Nt,Num_User,Num_BS))/sqrt(2);
         % New Small Scale Fading
         H = t_corr.*H + sqrt(1-t_corr^2).*N;
end