% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/08

%% function used to generate the optimal precoder by using KKT optimsiation form

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - h_k               (NtxN_UE double) = The effective channel for users
% - SINR              (1x1 double) = QoS
% - sigma             (1x1 double) = Noise power
% - tolerance         (1x1 double) = The accuracy that require for the
% lambda
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output
% - power (N_UEx1 double) = Transmit power for the transmitter to satisfy
% the QoS for all users.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Reference: Wireless communication and optimisations lecture note Part4(optimisation, KKT)
% 2024
function [power]=KKT_optimisation_RIS(h_k,SINR,sigma,tolerance)
        [Nt,N_UE] = size(h_k);
        % Fixed-point algorithm -> find the optimal lambda for the user
        % channel

        % initilise the lambda_old as 0
        lambda_old = zeros(N_UE,1);
        % the initile lambda input is 0.
        lambda = zeros(N_UE,1);
        while(1)
            % assign each loop the all user power matrix is zero
            all_ue_power = 0;
            for i = 1:N_UE 

                % generate the channel hermitiam for each user
                all_ue_power = all_ue_power + (lambda(i)/sigma)*h_k(:,i)*h_k(:,i)';
            end
    
            for i = 1:N_UE
                % user channel
                H_UE = h_k(:,i);

                % compute the user interference hermitian channel
                interference = all_ue_power - (lambda(i)/sigma)*(H_UE*H_UE');

                % compute the new lambda
                lambda(i) = real(SINR*sigma *( H_UE'*(eye(Nt) + interference )^(-1) *H_UE )^(-1));
            end

            % compute the difference between the old and new lambda, to
            % make sure the lambda value is converged to a optimal lambda
            % value
            accuracy = sum(lambda_old - lambda);

            % jump out the loop when the accuracy is less than the
            % tolerance
            if (abs(accuracy) < tolerance)
                break;
            end

            % update the old lambda to new lambda if the tolerance is not
            % statisfy
            lambda_old = lambda;
        end
        

        % Compute the optimal precoder based on the optimal lambda

        % compute the all user power
        ue_optimal_power = 0;
        for i = 1:N_UE
            H_UE = h_k(:,i);
            ue_optimal_power =  ue_optimal_power + (lambda(i)/sigma).*(H_UE*H_UE');
        end

        w_optimal = zeros(Nt,N_UE);
        for i = 1:N_UE
            H_UE = h_k(:,i);
            % Find the interference power
            interference = ue_optimal_power - (lambda(i)/sigma).*(H_UE*H_UE');
            % Compute the F for the precoder
            F = (eye(N_UE) + interference )^(-1) * H_UE;

            % Precoder direction
            w_optimal(:,i) = F/norm(F);
        end

        % Find the transmit power by forming the M matrix
        A = abs(h_k'*w_optimal).^2;
        B = diag(diag(A)) - A;
        M_matrix = B + (1/SINR).*(diag(diag(A)));
        M_inverse = (M_matrix)^(-1);
        power = M_inverse * (sigma.*ones(N_UE,1));
end