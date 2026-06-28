%% SDP optimisation for minimisation transmit power 
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
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [Sum_B_trace]=SDP_optimisation(H,sigma,SINR)

         % Find the parameters
         [Nr,Nt,Num_UEs] = size(H);
         H_new = squeeze(H); % H_new(:,1) channel for user 1

         HH = zeros(Nt,Nt,Num_UEs);
         % find the channel matrix h*h^(H) which should a matrix have
         % size(Nt*Nt)
         for i=1:Num_UEs
             HH(:,:,i) = H_new(:,i)*H_new(:,i)';
         end

         % Perfrom SDP optimisation
         cvx_begin quiet
         variable B(Nt,Nt,Num_UEs) hermitian % the precoder hermititan matrix B=w*w^(H), have size Nt*Nt*Num_UEs
         variable s(Num_UEs,1)  % the value that make the inequality become equality

         % Find trace for each hermitian function
         Sum_B_trace=0;
         for i = 1:Num_UEs
             % Accumulate the trace of the B hermitian matrix
             Sum_B_trace = Sum_B_trace+trace(B(:,:,i));
         end

         % Objective
         minimise Sum_B_trace
         % constraint
         subject to
               for i = 1:Num_UEs
    
                   % Target user channel trace
                   UE_trace = trace(HH(:,:,i)*B(:,:,i));
    
                   sum_BQ_ij = 0;
                   for j = 1:Num_UEs
                      if (i~=j)
                         % Compute other precoder matrix with the target user channel
                         sum_BQ_ij = sum_BQ_ij + trace(HH(:,:,i)*B(:,:,j));
                      end
                   end
    
                   % First constraint - s to become a equality constraint
                   UE_trace - SINR*sum_BQ_ij - s(i,1) == SINR*sigma
                   % Make sure is >=0, such that the first constraint is
                   % equality constraint  1)A>=0 2)A-S=0, So 3)S>=0
                   s(i,1)>=0
    
                   % Make sure the B matrix (w*w^(H)) is semidenfinitie, to
                   % follow the SDP optimisation
                   B(:,:,i) == hermitian_semidefinite(Nt)              
               end
         cvx_end
end