clc
clear
close all;
addpath("function_Part2\");
addpath("function_Part1\")
orange = [0.8500 0.3250 0.0980];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% Rate Maximisation 

%% system parameters
sigma = 1;                        % The noise power for each user is the same
Nt = 4;                           % The number of transmit antenna
Nr = 1;                           % Each user have single receive antenna
SNR_dB = -10:5:30;                % The snr range in dB
SNR_all = 10.^(SNR_dB./10);       % The snr value in linear scale
Channel_max = 100;                % The Number of channel realisation 
Num_UEs = 4;                      % The Number of user in the system
p_tolerance = 1e-10;              % The tolerance value for the water-filling
LW = 1.8;                         % The line width for plotting

C_sum_ZFWF = zeros(Channel_max,length(SNR_all));
C_sum_MRT  = zeros(Channel_max,length(SNR_all));
C_sum_MMSE = zeros(Channel_max,length(SNR_all));

i_max = 0;
while(i_max < Channel_max)
   i_max = i_max + 1;
  
   % Generate the channel for all users
   H = (randn(Nr,Nt,Num_UEs) +1i*randn(Nr,Nt,Num_UEs))./sqrt(2);
     for i_snr = 1:length(SNR_all)
        SNR = SNR_all(i_snr);

        
        % H(:,:,1) for user 1 with size Nr*Nt
% 
%         H_new = squeeze(H).';
%         F = H_new'*(H_new*H_new')^(-1);
%         W_ZFBF = zeros(Nt,Nt);
%         for i = 1:Num_UEs
%             W_ZFBF(:,i) = F(:,i)/norm(F(:,i));
%         end
%         User_power = H_new*W_ZFBF;
% 
%         [~,S,V]=svd(User_power);
%         lambda = diag(S.^2);
% 
%         % Iterative water filling algorithm (power allocation)
%         i = 0;                             % initialise the count
%         % loop tp find the constant value miu, that allocate power to all
%         % channels with non-negative power.
%         while(true)
%             i = i+1;                        % +1 for each iteration
%             inv1 = 1./(SNR.*lambda);        % compute the inverse of channel and SNR
% 
%             % compute the water level
%             u = (1/(Num_UEs-i+1))*(1+sum(inv1(1:Num_UEs-i+1)));
% 
%             % allocate power to channels
%             power = max(u - inv1,0);
% 
%             % The 'if' statement is used to determine whether the allocated power is
%             % under power constarint, usually it 1. (whether it is converge)
%             if ( abs(sum(power)-1) < p_tolerance)
%                     break;
%             end
%         end

        [power,lambda]=ZF_WF(H,SNR,p_tolerance);
        C = log2(1+(SNR.*power).*lambda);



        H_new = squeeze(H).';
        % equal power allocation
        P_equal = sqrt((1/Num_UEs)*(SNR*sigma));

        %% MMSE precoder (normalised)
        
%         W_MMSE = zeros(Nt,Num_UEs);
%         C_MMSE = zeros(1,Num_UEs);
% 
%         F_MMSE = H_new'*(SNR*(H_new*H_new') + eye(Num_UEs))^(-1);
%         for i=1:Num_UEs
%             W_MMSE(:,i) = P_equal*(F_MMSE(:,i)/norm(F_MMSE(:,i)));
%         end
%         A = abs(H_new*W_MMSE).^2;
%         for i = 1:Num_UEs
%             UE_p = A(i,i);
%             inter_MMSE = sum(A(i,:)) - UE_p;
%             C_MMSE(i) = log2(1+UE_p/(inter_MMSE+sigma));
%         end
        [C_MMSE]=W_MMSE(H,SNR,sigma,P_equal);


        
        %% MRT
        W_MRT = zeros(Nt,Num_UEs);
        C_MRT = zeros(1,Num_UEs);
        for i =1:Num_UEs
            h_ue = H(:,:,i);
            W_MRT(:,i) = P_equal*(h_ue'/norm(h_ue));
        end
        B = abs(H_new*W_MRT).^2;

        for i = 1:Num_UEs
            UE_power_MRT = B(i,i);
            MRT_inter = sum(B(i,:)) - UE_power_MRT;
            C_MRT(i) = log2(1+UE_power_MRT/(MRT_inter+sigma));
        end
        
        C_sum_ZFWF(i_max,i_snr) = sum(C);
        C_sum_MMSE(i_max,i_snr) = sum(C_MMSE);
        C_sum_MRT(i_max,i_snr) = sum(C_MRT);

     end
end
Avg_C_ZFWF = mean(C_sum_ZFWF,1);
Avg_C_MMSE = mean(C_sum_MMSE,1);
Avg_C_MRT  = mean(C_sum_MRT,1);

% Multiplexing gain
ZF_WF_multiplexing = Avg_C_ZFWF./log2(SNR_all)
MMSE_multiplexing = Avg_C_MMSE./log2(SNR_all)
MRT_multiplexing = Avg_C_MRT./log2(SNR_all)
figure()
plot(SNR_dB,Avg_C_ZFWF,"r-x",'LineWidth',LW,'MarkerSize',7);
hold on;grid on;
plot(SNR_dB,Avg_C_MMSE,'b-o','LineWidth',LW,'MarkerSize',7);
hold on;
plot(SNR_dB,Avg_C_MRT,'k-V','LineWidth',LW,'MarkerSize',7);
xlabel("SNR [dB]");ylabel("Average sum-rate [bps/Hz]")
legend("ZF-WF","MMSE (Regularised ZF-BF)","MRT",'Location','northwest');
title("Rate Maximisation by ZF-WF,MMSE,MRT");

figure()
plot(SNR_all,Avg_C_ZFWF,"r-x",'LineWidth',LW,'MarkerSize',7);
hold on;grid on;
plot(SNR_all,Avg_C_MMSE,'b-o','LineWidth',LW,'MarkerSize',7);
hold on;
plot(SNR_all,Avg_C_MRT,'k-V','LineWidth',LW,'MarkerSize',7);
xlabel("SNR [linear]");ylabel("Average sum-rate [bps/Hz]")
legend("ZF-WF","MMSE (Regularised ZF-BF)","MRT",'Location','northwest');
title("Rate Maximisation by ZF-WF,MMSE,MRT");
