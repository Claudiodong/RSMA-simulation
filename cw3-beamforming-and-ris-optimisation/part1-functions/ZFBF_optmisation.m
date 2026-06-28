%% ZFBF optimisation for minimisation transmit power 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Input
% - H (Nr*Nt*Num_UEs, complex) = User channel
%
% - sigma (1x1, double) = noise power for the user, consider as 1 for all
% in the simulation
%
% - SINR (1x1, double) = Required SINR value for the users (Quality of service)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output
% - tau (1x1, double) = Required Transmit power at Base Station for all
% users
%
% In the function, the ZFBF transmit precoder is designed to beamforming
% the direction, and the power is designed based on the precoder direction 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [tau]=ZFBF_optmisation(H,SINR,sigma)
       [Nr,Nt,Num_UEs] = size(H);
       % ZFBF precoder design
       H_new = squeeze(H).'; 
       F = H_new'*(H_new*H_new')^(-1);
       W_ZFBF = zeros(Nt,Num_UEs);
       for i=1:Num_UEs
       %    F = H(:,:,i)'*(H(:,:,i)*H(:,:,i)')^(-1);
           W_ZFBF(:,i) = F(:,i)/norm(F(:,i));
%           W_ZFBF(:,i) = F/norm(F);
       end
       

       cvx_begin quiet

           % Define the power variable, which is always positive number
           variable power_ZF(Num_UEs,1) nonnegative
           variable tau nonnegative

           % objective, minimise the transmit power
           minimise tau

           subject to
           for i=1:Num_UEs
               sum_users = 0;
               UE_H = H_new(i,:);
               for j = 1:Num_UEs
                   if (i~=j)
                    sum_users = sum_users + power_ZF(j,1)*square_abs(UE_H*W_ZFBF(:,j));
                   end
               end
                (sum_users+sigma) <= (1/SINR)*(power_ZF(i,1)*square_abs(UE_H*W_ZFBF(:,i)))
           end
           sum(power_ZF) <= tau

       cvx_end

end