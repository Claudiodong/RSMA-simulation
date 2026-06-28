clc
clear
close all;
addpath("New Function\");
orange = [0.8500 0.3250 0.0980];
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex'); % For text
set(0, 'defaultLegendInterpreter', 'latex'); % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((N0_dBm-30)/10);               % noise variance in linear
Num_User = 10;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
Nt1 = [16,64];                          % Number of Tx antenna
%Nt1 = [4,16,64];                        % Number of Tx antenna
Nr = 2;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels (time correlation)
time_k = 10;                           % time instant (spatial correlation)
drop = 1e3;                             % Number of drop
Q = 4;                                  % the number of scheduled user at this time instant
Avg_C_time = zeros(drop,Num_User);
Avg_C_sort = zeros(length(Nt1),drop*Num_User);

for i_nt = 1:length(Nt1)
    Nt = Nt1(i_nt);
    % drops loop, where the UEs locations unchange,
    for i_drop = 1:drop
    
         % center BS to users
         [User_location,center_BS] = Random_user_deploy(d(2),d(1),Num_User,center);
                
         % interfer BS location
         BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
         BS_location = zeros(2,Num_BS+1);
         Distance = zeros(Num_BS+1,Num_User);
    
         for i = 2:Num_BS+1
              % interfere BSs location 
              BS_location(:,i) = r*[cos(BS_theta(i-1));sin(BS_theta(i-1))];
    
              % Distance 
              Distance(i,:) = sqrt(sum((BS_location(:,i) - User_location).^2));
         end
         Distance(1,:) = center_BS;
    
         % Compute the path loss in dB
         PL_dB = 128.1+37.6*log10(Distance./1000)+shadowing_dB*randn(Num_BS+1,Num_User);
         % Path loss in linear, Center BS to UEs
         PL = 10.^((PL_dB)./10);
         
         % initilise the small scale channel 
         %[~,H] = small_scale_channel(Nr,Nt,Num_User,Num_BS+1,time_k,t_corr);
         % Time instant loop where the small scale fading change
         H = (randn(Nr,Nt,Num_User,Num_BS+1) + 1i* randn(Nr,Nt,Num_User,Num_BS+1))/sqrt(2);
    
    
         % Uniform Power Allocation 
         Ne = min(min(Nt,Num_User),Q);
         Spower = Tx_power/Q/Nr;
         
         Sum_C = zeros(time_k,Num_User);
         % time instant loop
         for i_t = 1:time_k
             % Random User Schedule
             
             % The index of the scheduled User
             Scheduled_User = randperm(Num_User,Q);
    
             % Scheduled User Channels for all BS(including the center BS)
             [Scheduled_H]=Scheduled_User_H(Scheduled_User,H,PL);
    
             % Separate the center and interfere channels
             center_H = (Scheduled_H(:,:,:,1)); 
             H_interfere = ((Scheduled_H(:,:,:,2:end)));
    
             % Random precoder design for the interfere BSs
             w_interfere = (randn(Nt,Nr,Q,Num_BS) + 1i*randn(Nt,Nr,Q,Num_BS))./sqrt(2) ;
             % normalise the precoder for the interfere BS
             w_interfere = w_interfere./pagenorm(w_interfere);
    
             % Define the size of the I_o
             I_o = zeros(Nr,Q,Num_BS);
             R_noise = zeros(Nr,Nr,Q);
    
             % find the inter-cell interference for covariance matrix
             for i = 1:Q
                 index = 2*i-1:2*i;
                 for j_bs = 1:Num_BS
                     UE_channel = H_interfere(:,:,i,j_bs);
                     w_UE = w_interfere(:,:,i,j_bs);
                     L = UE_channel*w_UE;
                     % interference from interference user
                     I_o(:,index,j_bs) = Spower.*(L*L');
                     % I_o = [H1*W1,H1*W2 , H2*W1,H2*W2] for every interfere BS
                 end
                 % compute the noise covaraince matrix 
                 I_o_UE =  sum(I_o(:,index,:),3);  
                 R_noise(:,:,i) = I_o_UE+N0*eye(Nr);
             end
                
%              for i = 1:Q
%                   index = 2*i-1:2*i;
%                   % inter-cell interference ( interfere BS from other cell to scheduled user)
%                   I_o_UE =  sum(I_o(:,index,:),3);  
%                   R_noise(:,:,i) = I_o_UE+N0*eye(Nr);        
%              end
            
             % Precoder design for the whitening channel
             [w,g]=Block_diagonalization_design(center_H,R_noise,Spower);
  
             A = zeros(Nr*Q,Nr*Q);
             for i = 1:Q
                 index = 2*i-1:2*i;           
                 for j = 1:Q  
                     G = g(:,:,i);
                     O = G*center_H(:,:,i)*w(:,:,j);
                     % Scheduled user power matrix (H1_1*W1_1,H1_1*W1_2;H1_2*W1_1,H1_2*W1_2)
                     %A(index,(2*j-1):2*j) = O*O';   
                     A(index,(2*j-1):2*j) = abs(O).^2;    
                 end    
             end
    
             % compute the SINR
             for i_q = 1:Q
                 index = 2*i_q-1:2*i_q;
                 % user power
                 User_Power = A(index,index);
                 % all the inter-stream interference
                 I_l = sum(User_Power.',2) - diag(User_Power);
    
                 % Intra-cell interference (co-scheduled user interference)
                 I_c = sum(A(index,:),2) - sum(User_Power,2);
    
                 % Combiner
                 G = g(:,:,i_q);

                 % inter-cell interference from interfere BSs
                 I_o1= 0;
                 for i_BS = 1:Num_BS
                     UE_channel1 = H_interfere(:,:,i_q,i_BS);
                     % Assign the precoder with uniform power allocation
                     w_UE1 = w_interfere(:,:,i_q,i_BS)*(sqrt(Spower)*eye(Nr));
                     I_o1 = I_o1 + (G*UE_channel1*w_UE1)*(G*UE_channel1*w_UE1)';
                 end

                 % All the interference for the user q
                 interference = I_l + I_c + sum(I_o1,2) + vecnorm(G,2,2).^2 *N0;
    
                 % SINR
                 SINR = diag(User_Power)./(interference);
                 % Capacity
                 C = log2(1+SINR);
                 
                 % allocate the user rate into corresponding index 
                 Sum_C(i_t,Scheduled_User(i_q)) = sum(C);
             end    
             % Update the Small scale fading
             [H]=Update_Small_Scale_Fading(H,t_corr);             
         end
         % Average User Rate per drop after time_k time
         Avg_C_time(i_drop,:) = sum(real(Sum_C),1)./time_k;
    end
    % the average user rate for different Nt number
    Avg_C_sort(i_nt,:) = sort(reshape(Avg_C_time,1,drop*Num_User),"ascend");
    sprintf("\n Nt = %1.0f is simulated",Nt1(i_nt))
end
%%
% compute the CDF index x-axis
CDF = (1:length(Avg_C_sort))./length(Avg_C_sort)*100;

%%
figure()
for i = 1:length(Nt1)
    plot(Avg_C_sort(i,:),CDF);
    hold on;   
end
grid on;xlabel("Average User Rate [bps/Hz]");ylabel("Cumulative Distribution Function [\%]")
legend("$Nt=4$","$Nt=16$","$Nt=64$")

%save("Nr2Q4","Avg_C_sort");
