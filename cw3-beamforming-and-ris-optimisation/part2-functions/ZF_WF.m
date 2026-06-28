
%% Zero-forcing Beamforming with water-filling algorithm (ZF-WF)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - H (Nr*Nt*Num_UEs complex) = User Channel, iid Rayleigh fading channel (normalised)
%
% - SNR (1x1 double) = Signal to Noise ratio
%
% - p_tolerance = power tolerance value, to make sure the sum of all power
% is under constaint of 1.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Power (Num_UEs*1) = the allocated power to each user channel which
% experience interference free, and also to maximise the sum rate by
% allocate more power to strong channel and vice versa for weak channel
%
% - lambda (Num_UEs*1, double) = The transmission channel power (in this case
%  which is from the product of ZFBF precoder and channel)
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [power,lambda]=ZF_WF(H,SNR,p_tolerance)
      
        % define the parameter
        [Nr,Nt,Num_UEs] = size(H);
        % squeeze the channel into size (Num_UEs*Nt), H_new(1,:)= UE1
        % channel
        H_new = squeeze(H).';

        %% Gnerating the ZFBF
        % Find the F matrix of ZFBF 
        F = H_new'*(H_new*H_new')^(-1);
        W_ZFBF = zeros(Nt,Nt);
        for i = 1:Num_UEs
            % Normalised the ZFBF precoder
            W_ZFBF(:,i) = F(:,i)/norm(F(:,i));
        end

        % Perfrom the User with precoder matrix H1*W(1:K) and so on, should
        % be a diagnal matrix
        User_power = H_new*W_ZFBF;

        % Performing the SVD, and find the singular value
        [~,S,~]=svd(User_power);
        % Find the lambda which is squared of the singular value
        lambda = diag(S.^2);

        %% Iterative water filling algorithm (power allocation)
        i = 0;                             % initialise the count
        % loop tp find the constant value miu, that allocate power to all
        % channels with non-negative power.
        while(true)
            i = i+1;                        % +1 for each iteration
            inv1 = 1./(SNR.*lambda);        % compute the inverse of channel and SNR

            % compute the water level
            u = (1/(Num_UEs-i+1))*(1+sum(inv1(1:Num_UEs-i+1)));

            % allocate power to channels
            power = max(u - inv1,0);

            % The 'if' statement is used to determine whether the allocated power is
            % under power constarint, usually it 1. (whether it is converge)
            if ( abs(sum(power)-1) < p_tolerance)
                    break;
            end
        end

end