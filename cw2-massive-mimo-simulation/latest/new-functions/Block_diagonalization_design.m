%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Input
% - center_H (Nr*Nt*Q complex)= The center BS channel to scheduled users 
% - R_noise(Nr*Nr*Q complex) = The whitening noise
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
% - w_center(Nt*Nr*Q complex) = The w ZFBF precoder design for the center BS to all
% scheduled user
% - g (Nr*Nr complex) = The combiner design
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Using the function of block diagonlization which require the rank of the
% interference space is less than Nt(number of transmit antenna), also,
% K*Nr, where K is the number of interference user (Q-1), which K*Nr<Nt is
% also need to statisfied, otherwise, there is not null space.

% The design force W_p lie in the null space of the H_q (interference space),
% Therefore, the product of H_q*W_p is 0, which is orthogonal to each
% other.
%
% Method: 1) Find the interference space of the intended scheduled user
%
%         2) Perfrorm the SVD (Singular Vector Decompsition) on H_q_null
%            where [U,S,V]=svd(H_q_null) H_q_null = the interference space
%            of H_q
%
%         3) Find the number of non-zero singular value in S, Null_index =
%            (Number of non-zero singular +1).
%            V_H = [V_q,V_q_dot]_H
%
%         4) The null space is the eignvector of V_q_dot = V(:,Null_index+1) 
%            which its coressponding value in S is 0. 
%
%         5) Then find the dominant eignvector H_eq_q = H_q*V_q_dot
%            [~,S_eq,V_eq]=svd(H_eq_q), where V_eq_H = [V_eq_q,V_eq_q_dot]_H
%            H_q = 用户信道

%            V_eq_q = V_eq_H(:,1:Nr); The number of stream that require to transmit
%
%         6) Then the precoder w_q can be design as
%            w_q = V_q_dot*V_eq_q
%
% All of this need H_q satisfy two conditions
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [w,g] = Block_diagonalization_design(center_H,R_noise,Spower)
            
             [Nr,Nt,Q]=size(center_H);

             % Reduce the size of matrix and change it into (Nr*Q)*Nt form
             H_all = zeros(Nr*Q,Nt); % All Scheduled User channels
             for i_q = 1:Q
                 star = (i_q*2)-1;
                 en = i_q*2;
                 H_all(star:en,:) = center_H(:,:,i_q);
             end
    
             % The size of the precoder
             w = zeros(Nt,Nr,Q);

             % Precoder design by using block diagonalization knowledge
             for i_q = 1:Q
                 % Reset interference space for every loop start
                 H_q = H_all;  
                 star = (i_q*2)-1;
                 en = i_q*2;
                 % eliminate the user channel, left with interference space
                 H_q(star:en,:) = [];
                 % Where the number of interference user should be K = Q-1

                 % Check the rank of interference space, whether satisfy
                 % The condition " rank(H_interference_space) < Nt "
                 % Also K*Nr < Nt, otherwis there is not null space
                 condition_1 =  rank(H_q) > Nt ;
                 condition_2 = Nr*(Q-1) > Nt;

                 if ( condition_1 == 1 || condition_2 == 1)

                    % having error when the condition is not satisfied
                    error(sprintf("Condition 1 or 2 is not satisfied, " + ...
                        "the null space of H_q is 0 for this set up, Nt=%1.0f,Q=%1.0f",Nt,Q));
                 end
    
                 % Perform the SVD for the interference space
                 [~,S_s,V_s]=svd(H_q);

                 % Find the null space index by find the number of non-zero singular values +1 
                 Null_index = length(diag(S_s))+1;

                 % The null space for the H interference space
                 V_q_space = V_s(:,Null_index:end);

                 % whitening noise
                 Whitening = R_noise(:,:,i_q)^(-1/2);

                 % The optimal solution which transmit through the nth
                 % domiant eignvectors, with whiteneing
                 [U_eq,~,V_eq_q]=svd(Whitening*center_H(:,:,i_q)*V_q_space);
% 
%                  % Find the number of non-zero singular value
%                  N_index = length(diag(S_eq_q));

                 % Compute the precoder design based on interference space
                 % and effective domiant eignvector
                 w(:,:,i_q) = (V_q_space*V_eq_q(:,1:Nr))*sqrt(Spower)*eye(Nr);

                 % combiner design
                 g(:,:,i_q) = U_eq'*(R_noise(:,:,i_q)^(-1/2));
             end
end