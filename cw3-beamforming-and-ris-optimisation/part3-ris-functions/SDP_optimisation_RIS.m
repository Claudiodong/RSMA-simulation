% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/08

%% function used to generate the optimal precoder by using SDP optimsiation form
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - h_k               (NtxN_UE double) = The effective channel for users
% - SINR              (1x1 double) = QoS
% - sigma             (1x1 double) = Noise power
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Ouput
% - sum_B_trace (1x1 double) = total transmit power for the transmitter to satisfy
% the QoS
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function [sum_B_trace]=SDP_optimisation_RIS(h_k,SINR,sigma)
         % find the paramters
        [Nt,N_UE] = size(h_k);

        % Perfrom SDP optimisation
        cvx_begin quiet% the begin of the cvx optimisation
          
            % Define the HH hermititan user channel HH(:,:,1) for user 1
            HH = zeros(Nt,Nt,N_UE);
            for i = 1:N_UE
                % 1e3 is the scaling factor which can increase the
                % resolution

                % user hermitian channel matrix
                HH(:,:,i) = (h_k(:,i)*h_k(:,i)');
            end

            % Define the variables
            variable B(Nt,Nt,N_UE) hermitian % the precoder hermitian matrix
            variable s(N_UE,1)
            
           
            % Define the objective function
            sum_B_trace = 0;
            for i = 1:N_UE
                % Accumualte the trace for each B matrix 
                sum_B_trace = sum_B_trace + trace(B(:,:,i));
            end

            % Objective function, minimise the transmit precoder power
            minimise (sum_B_trace)

            % define the constraint
            subject to
            for i = 1:N_UE
                UE_power = trace(HH(:,:,i)*B(:,:,i));
                interference = 0;
                for j = 1:N_UE
                    if (i~=j)
                        interference = interference + trace(HH(:,:,i)*B(:,:,j));
                    end
                end

                % first constraint => make it as equality constraint
                % SINCE the HH is scaling up by 1e3, then the RHS should be
                % also scaled up by 1e3
%                 1e5*(UE_power - SINR*interference - s(i,1)) == 1e5*SINR*(sigma)
%                 
%                 % second constraint which make the first constraint as
%                 % equality constrain
%                 s(i,1) >= 0


                1e10*(UE_power - SINR*interference ) >= 1e10*SINR*sigma

                % Third constraint
                B(:,:,i) == hermitian_semidefinite(Nt) % The semitdefinite constraint
            end
        cvx_end % end of cvx optimisation

end