% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/06
% 2024/03/08 updated

%% function used to generate the optimal phase shift matrix from RIS to users
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Input
% - h_d_k              (NtxN_UE complex) = The direct channel from
% transmitter to users
% - h_r_k            (MxN_UEs complex) = The reflect channel from RIS elements to
% users
% - G              (MxNt complex) = The channel from transmitter to RIS
% elements
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output
% - O_optimal (M*M double) = The optimal phase shift from the RIS to
% users
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Reference:
% 1) Q. Wu and R. Zhang, 
% "Intelligent Reflecting Surface Enhanced Wireless Network: Joint Active and Passive Beamforming Design," 
% 2018 IEEE Global Communications Conference (GLOBECOM), Abu Dhabi, United Arab Emirates, 2018, pp. 1-6,
% doi: 10.1109/GLOCOM.2018.8647620.

% 2) Q. Wu and R. Zhang, 
% "Intelligent Reflecting Surface Enhanced Wireless Network via Joint Active and Passive Beamforming," 
% in IEEE Transactions on Wireless Communications, vol. 18, no. 11, pp. 5394-5409, Nov. 2019, 
% doi: 10.1109/TWC.2019.2936025.

function [O_optimal]=RIS_phase_optimisation(h_d_k,h_r_k,G)
       [M,~] =size(G);
       [Nt,N_UE]=size(h_d_k);
       %% Maximise the sum rate
       cvx_begin quiet

          variable V(M+1,M+1) hermitian
          R = zeros(M+1,M+1);
          sum_trace = 0;
          for i = 1:N_UE
              phi = (diag(h_r_k(:,i)')*G); %  diag(diag(h_r_k')) make the matrix N*N
              % phi should have size M*N
    
              % Sub the value into the matrix
              R(1:M,1:M) = phi*phi';
              R(1:M,M+1) = phi*h_d_k(:,i);
              R(M+1,1:M) = h_d_k(:,i)'*phi';
    
              sum_trace = sum_trace+real(trace(R*V)) + square_pos(norm(h_d_k(:,i)')) ;
          end
    
           maximise sum_trace
    
           subject to 
            % constraint
           
               for j = 1:M+1
                   (V(j,j)) == 1 % The diagonal should have modulo of 1
               end
               % The V matrix should be semidefinite
               V == hermitian_semidefinite(M+1)   
       cvx_end

       %% using eign value decompsitoin to find the rank 1 solution
%        O_optimal = zeros(M,M,N_UE);
%        for i = 1:N_UE
%            [U,lambda]=eig(V(:,:,i));
%            % Find the largest eignvalue index
%            r = (randn(M+1,1) +1i* randn(M+1,1)) ./sqrt(2);
%            v_hat = U*sqrt(lambda)*r;
%            % convert the angle to radians
%            v_hat_M = angle(v_hat(1:M)./v_hat(M+1));
%            % convert into complex format
%            v = exp(1i*v_hat_M);
%            % make it as diagonal matrix
%            O_optimal(:,:,i) = diag(v);
%        end
       
       % Do the eignvalue decomposition
       [U,lambda] = eig(V); % the solution is usually rank 2, not rank 1 by checking the lambda value
       % find the maximum eignvalue index
       [~,index] = max(diag(lambda));
       % find the eignvector
       v_hat = U(:,index);
       % find the angle, and need to hermitian transpose
       da_V = exp(1i*angle(v_hat(1:M)'));
       O_optimal = diag(da_V);       
end