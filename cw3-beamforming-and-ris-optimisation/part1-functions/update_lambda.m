
%% Update the Lagrangian Multipler lambda by using fixed-point Method
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - H_all (Num_UEs*Nt,complex) = User channel H(:,1) = user channel for
% user 1
% - Lambda (Num_UEs*1, double) = lagrangian multiplex, which must be >=0,
% such that to reflect the sensitive of the channel.
% - sigma (1x1, double) = Noise power at user, in this simulation, the noise power
% is 1
% - SINR (1x1, double) = The signal-to-interference-plus-noise ratio
%
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% - Lambda_new (Num_UEs*1,double) = New lambda value
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function [Lambda_new] = update_lambda(H_all,Lambda,sigma,SINR)

            [Num_UEs,Nt]=size(H_all);

            HH_all= zeros(Nt,Nt);
            for j = 1:Num_UEs
                % the user channel
                H_ue = H_all(:,j);
                % compute the sum of all channels with lambda value
                HH_all = HH_all+Lambda(j)*(H_ue*H_ue');
            end
            
            Lambda_new = zeros(Num_UEs,1);

            % Update the lambda
            for i = 1:Num_UEs
                % user channel
                H_ue = H_all(:,i);
                % substract the target user from the sum channel matrix
                HH_sum = HH_all - Lambda(i)*(H_ue*H_ue');
                % calculate the new lambda based on previous lambda value
                Lambda_new(i) = real((sigma.*SINR)/(H_ue'*inv(eye(Num_UEs)+HH_sum)*H_ue));
                
            end

end