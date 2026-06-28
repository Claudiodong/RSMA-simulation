clc
clear
close all;
addpath("function_Part1\")
orange = [0.8500 0.3250 0.0980];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels

%% system parameters
sigma = 1;                        % The noise power for each user is the same
Nt = 4;                           % The number of transmit antenna
Nr = 1;                           % Each user have single receive antenna
SINR_dB = -5:3:10;                % The sinr range in dB
SINR_all = 10.^(SINR_dB./10);     % The sinr value in linear scale
Channel_max = 1;                % The Number of channel realisation 
Num_UEs = 4;                      % The Number of user in the system
tolerance = 1e-10;                % Tolerance for lagrangian multiplier convergen

Sum_Power_KKT  = zeros(Channel_max,length(SINR_all));
Sum_Power_SOCP = zeros(Channel_max,length(SINR_all));
Sum_Power_SDP  = zeros(Channel_max,length(SINR_all));
Sum_Power_ZFBF = zeros(Channel_max,length(SINR_all));


i_max = 0;% initilise the counter
% channel realisation while loop
while (i_max < Channel_max)
        i_max = i_max +1; % counter increment

    
        % Generate the channel for all users
        H = (randn(Nr,Nt,Num_UEs) +1i*randn(Nr,Nt,Num_UEs))./sqrt(2);
        % H(:,:,1) for user 1 with size Nr*Nt

    for i_snr = 1:length(SINR_all)
        SINR = SINR_all(i_snr); % The SINR value 
%         % perfrom SOCP optimisation with CVX toolbox help
%         H_new = squeeze(H).'; % H(1,:) = channel for user 1
% 
%         cvx_begin quiet
% 
%             % Define thr variable for the problem
%             variable tau nonnegative  % To make sure the transmit power is always positive
%             variable W(Nt,Num_UEs) complex % Define the precoder size with Nt*Num_UEs
% 
%             minimize(tau)
%             subject to
%             for i=1:Num_UEs
%                 norm( H_new(i,:) * W + sigma ) <= sqrt(1 + (1/SINR)) * real(H_new(i,:)*W(:,i))
%             end
%             norm(vec(W)) <= sqrt(tau)
% 
%         cvx_end

        % KKT optimisation 
        [Power_KKT]=KKT_optimisation(H,tolerance,sigma,SINR);
        % SOCP optimisation 
        [Power_SOCP,W]=SOCP_optimisation(H,sigma,SINR);
        % SDP optimisation 
        [Power_SDP]=SDP_optimisation(H,sigma,SINR);

        % Accumulating the power
        Sum_Power_KKT(i_max,i_snr)  = sum(Power_KKT);
        Sum_Power_SOCP(i_max,i_snr) = Power_SOCP;
        Sum_Power_SDP(i_max,i_snr)  =  Power_SDP;

%          % Perfrom SDP optimisation
%          cvx_begin quiet
%          variable B(Nt,Nt,Num_UEs) hermitian % the precoder hermititan matrix B=w*w^(H), have size Nt*Nt*Num_UEs
%          variable s(Num_UEs,1)
%          % find trace for each hermitian function
%          Sum_B_trace=0;
%          for i = 1:Num_UEs
%              % accumulate the trace of the B hermitian matrix
%              Sum_B_trace = Sum_B_trace+trace(B(:,:,i));
%          end
% 
%          % objective
%          minimise Sum_B_trace
%          % constraint
%          subject to
%            for i = 1:Num_UEs
%                % Target user channel trace
%                UE_trace = trace(HH(:,:,i)*B(:,:,i));
%                sum_BQ_ij = 0;
%                for j = 1:Num_UEs
%                   if (i~=j)
%                      % compute other precoder matrix with the target user channel
%                      sum_BQ_ij = sum_BQ_ij + trace(HH(:,:,i)*B(:,:,j));
%                   end
%                end
%                % first constraint - s to become a equality constraint
%                UE_trace - SINR*sum_BQ_ij - s(i,1) == SINR*sigma
%                % make sure is >=0, such that the first constraint is
%                % equality constraint  1)A>=0 2)A-S=0, So 3)S>=0
%                s(i,1)>=0
%                B(:,:,i) == hermitian_semidefinite(Nt)              
%            end
% 
% 
%          cvx_end
          
      [tau1]=ZFBF_optmisation(H,SINR,sigma);

       Sum_Power_ZFBF(i_max,i_snr)  = tau1;


%         % ZFBF precoder design
%        H_new = squeeze(H).'; 
%        F = H_new'*(H_new*H_new')^(-1);
%        W_ZFBF = zeros(Nt,Num_UEs);
%        for i=1:Num_UEs
%            W_ZFBF(:,i) = F(:,i)/norm(F(:,i));
%        end
% 
% 
%        cvx_begin quiet
%        variable p_ue(Num_UEs) nonnegative
%        variable tau nonnegative
%        
%        minimise (tau)
% 
%            UE_power = cvx(zeros(Num_UEs,Num_UEs));
%            for i=1:Num_UEs
%                for j = 1:Num_UEs
%                 UE_power(i,j) = p_ue(i)*square_abs(H_new(i,:)*W_ZFBF(:,j));
%                end
%            end
%     
%            subject to
%            for i=1:Num_UEs
%                sum_other_ue_p = sum(UE_power(i,:));
%                sum_other_ue_p + sigma <= (1+1/(SINR))*UE_power(i,i);         
%            end
%            sum(p_ue) <= tau
%            
% 
%        cvx_end
% 
%        sum_tau = sum_tau +tau;
    end
   
end
% Average transmit power
Avg_power_KKT  = mean(Sum_Power_KKT,1);
Avg_power_SOCP = mean(Sum_Power_SOCP,1);
Avg_power_SDP  = mean(Sum_Power_SDP,1);
Avg_power_ZFBF = mean(Sum_Power_ZFBF,1);
% convert from linear to dB scale
Avg_power_dB_KKT = 10.*log10(Avg_power_KKT);
Avg_power_dB_SOCP = 10.*log10(Avg_power_SOCP);
Avg_power_dB_SDP = 10.*log10(Avg_power_SDP);
Avg_power_dB_ZFBF = 10.*log10(Avg_power_ZFBF);

%%
figure()
%subplot(2,1,1)
plot(SINR_dB,Avg_power_dB_KKT,"b-*",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_dB_SOCP,"r-O",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_dB_SDP,"-X","Color",orange,'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_dB_ZFBF,"k-X",'LineWidth',1.5,'MarkerSize',10);
grid on;xlabel("SINR requirement [dB]");ylabel("Average Transmit Power [dB]");
%title("optimisation");
legend("KKT","SOCP","SDP","ZFBF",'Location','southeast')

figure()
% subplot(2,1,2)
plot(SINR_dB,Avg_power_KKT,"b-*",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_SOCP,"r-O",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_SDP,"-X","Color",orange,'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(SINR_dB,Avg_power_ZFBF,"k-X",'LineWidth',1.5,'MarkerSize',10);
grid on;xlabel("SINR requirement [dB]");ylabel("Average Transmit Power [W]");
%title(" optimisation");
legend("KKT","SOCP","SDP","ZFBF",'Location','northwest')

figure()
plot(1:Channel_max,10.*log10(Sum_Power_SOCP(:,1)),"g-*",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(1:Channel_max,10.*log10(Sum_Power_SDP(:,1)),"r-O",'LineWidth',1.5,'MarkerSize',10);
hold on;
plot(1:Channel_max,10.*log10(Sum_Power_KKT(:,1)),"-X","Color",orange,'LineWidth',1.5,'MarkerSize',10);
grid on;xlabel("Number of Channel Realisation");ylabel("Average Transmit Power [dB]");
%title("optimisation");
legend("SOCP","SDP","KKT",'Location','southeast')
