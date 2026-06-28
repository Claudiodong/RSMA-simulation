
%% SOCP optimisation for minimisation transmit power 
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
% - W (Nt*Num_UEs, complex) = transmit precoder design for the system
% Used to check the performance between the ZFBF and the optimisation
% methods.
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [tau,W]=SOCP_optimisation(H,sigma,SINR)

        [Nr,Nt,Num_UEs]=size(H);
        % perfrom SOCP optimisation with CVX toolbox help
        H_new = squeeze(H).'; % H(1,:) = channel for user 1

        cvx_begin quiet

            % Define thr variable for the problem
            variable tau nonnegative  % To make sure the transmit power is always positive
            variable W(Nt,Num_UEs) complex % Define the precoder size with Nt*Num_UEs

            minimize(tau)
            subject to
            for i=1:Num_UEs
                norm( [H_new(i,:) * W  sigma] ) <= sqrt(1 + (1/SINR)) * real(H_new(i,:)*W(:,i))
                imag(H_new(i,:)*W(:,i))==0
            end
            norm(vec(W)) <= sqrt(tau)

        cvx_end
end