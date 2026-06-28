% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/06

%% function used to generate the optimal precoder by using SOCP optimsiation form
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Input
% - h_k               (NtxN_UE double) = The effective channel for users
% - SINR              (1x1 double) = QoS
% - sigma             (1x1 double) = Noise power
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output
% - tau (1x1 double) = total transmit power for the transmitter to satisfy
% the QoS
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function [tau]=SOCP_optimisation_RIS(H,SINR,sigma)
       
        [Nt,Num_UEs]=size(H);

        cvx_begin quiet

            % Define thr variable for the problem
            variable tau nonnegative  % To make sure the transmit power is always positive
            variable W(Nt,Num_UEs) complex % Define the precoder size with Nt*Num_UEs

            minimize(tau)
            subject to
            for i=1:Num_UEs
                1e10*norm( [H(:,i)' * W  sqrt(sigma)] ) <= 1e10*sqrt(1 + (1/SINR)) * real(H(:,i)'*W(:,i))
                imag(H(:,i)'*W(:,i))==0
            end
            norm(vec(W)) <= sqrt(tau)

        cvx_end
end