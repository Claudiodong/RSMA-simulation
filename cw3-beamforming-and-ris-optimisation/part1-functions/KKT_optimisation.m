%% KKT optimisation for minimise the transmit power
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - H(Nr*Nt*Num_UEs,complex) = User channel H(:,:,1) = user channel for
% user 1
% - tolerance (1*1, double) = The tolerance value for the converge.
% - sigma (1x1, double) = Noise power at user, in this simulation, the noise power
% is 1
% - SINR (1x1, double) = The signal-to-interference-plus-noise ratio
%
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% - Power (Num_UEs*1,double) = The minimum transmit power required for
% achieving the target SINR value.
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [Power]=KKT_optimisation(H,tolerance,sigma,SINR)
        % find the parameters
        [Nr,Nt,Num_UEs]=size(H);

        H_all = squeeze(H);% squeeze the size into Num_User*Nt

        % initilise the lambda value
        Lambda_old = 0; % the old lambda for accuracy calacualtion
        Lambda_new = zeros(Num_UEs,1);
        
        % Fixed-point iteration
        while(1)
%             HH_all= zeros(Nt,Nt);
%             for j = 1:Num_UEs
%                 H_ue = H_all(:,j);
%                 HH_all = HH_all+Lambda_new(j)*(H_ue*H_ue');
%             end
%             
%             % update the lambda
%             for i = 1:Num_UEs
%                 H_ue = H_all(:,i);
%                 HH_sum = HH_all - Lambda_new(i)*(H_ue*H_ue');
%                 Lambda_new(i) = real((sigma*SINR)/(H_ue'*inv(eye(Num_UEs)+HH_sum)*H_ue));
%             end

            [Lambda_new] = update_lambda(H_all,Lambda_new,sigma,SINR);

            % The accuracy, make sure it is converged
            accuracy = sum(abs(Lambda_new-Lambda_old));

            % Break when below the tolerance
            if (accuracy<tolerance)
                 break;
            end

            % Update the lambda_old
            Lambda_old = Lambda_new;
        end
        
        HH_all= zeros(Nt,Nt);
        for j = 1:Num_UEs
            % perfrom the new lambda_value with channel matrix
            HH_all = HH_all+Lambda_new(j)*(H_all(:,j)*H_all(:,j)');
        end

        % compute the precoder design
        W = zeros(Nt,Nt);
        for j = 1:Num_UEs
            H_ue = H_all(:,j);
            HH_sum = HH_all - Lambda_new(j)*(H_ue*H_ue');
            % compute the F matrix
            F = inv(eye(Nt)+HH_sum) * H_ue;
            % Normalised the precoder
            W(:,j) = F/norm(F);
        end

        % All channel power
        HH_power = abs(H_all'*W).^2;
        % substract the diagonal power,and negative the other such as
        % h(j)w(k)
        B = diag(diag(HH_power)) - HH_power;
        % compute the M channel
        M = (1/SINR)*diag(diag(HH_power)) + B;
        % The inverse M matrix
        M_inv = M^(-1);
        % Compute the Transmit power based on the M
        Power = M_inv*ones(Num_UEs,1);

        % Check user transmit power, whether have error such that negative
        % power which is impossible
        if(Power<0)
            error("Negative power")
        end

end