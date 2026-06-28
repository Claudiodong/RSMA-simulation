function [C_MMSE]=W_MMSE(H,SNR,sigma,P_equal)

        % Find the parameter
        [Nr,Nt,Num_UEs] = size(H);

        % squeeze the channel since the Nr=1
        H_new = squeeze(H).';

        % define the variables
        W_MMSE = zeros(Nt,Num_UEs);
        C_MMSE = zeros(1,Num_UEs);

        % The F matrix for the MMSE
        F_MMSE = H_new'*(SNR*(H_new*H_new') + eye(Num_UEs))^(-1);
        % Normalised the MMSE precoder, and allocate the power on it 
        % (Euqal power allocation)
        for i=1:Num_UEs
            W_MMSE(:,i) = P_equal*(F_MMSE(:,i)/norm(F_MMSE(:,i)));
        end

        % Perfrom the User with precoder matrix H1*W(1:K) and so on, should
        % be a diagnal matrix
        A = abs(H_new*W_MMSE).^2;

        % Calculate the capacity
        for i = 1:Num_UEs
            UE_p = A(i,i);
            inter_MMSE = sum(A(i,:)) - UE_p;
            C_MMSE(i) = log2(1+UE_p/(inter_MMSE+sigma));
        end

end