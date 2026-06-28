clc
clear
close all;
orange = [0.8500 0.3250 0.0980];
%% parameter
Tx_power_dBm = 46;                      % Transmit power in dBm
Tx_power = 10^((Tx_power_dBm-30)/10);   % Transmit power in linear
N0_dBm = -174;                          % Noise variance in dBm
N0 = 10^((-174-30)/10);                 % noise variance in linear
Num_User = 10;                          % Number of user in the system
shadowing_dB = 8;                       % Log-normal shadowing in dB
shadowing_correlation = 0;
Nt = 4;                                 % Number of Tx antenna
Nr = 1;                                 % Number of receive antenna per user
d = [35 250];                           % minimum and maximum distance away from the BS
r = 500;                                % radius of the circle
Num_BS = 6;                             % number of other BS as interfer
center = [0,0];                         % center BS location
t_corr = 0.85;                          % time correlation for channels
time_k = 2;                             % time instant
MAX = 1;
Q = 4;                                  % random user s
i_runs = 0;
    while (i_runs<MAX)
        i_runs = i_runs +1; % counter for the monte carlo
    
        %% Random deploy User
        % center BS to users
        [User_location,Center_BS_UEs] = Random_user_deploy(d(2),d(1),Num_User,center);
        
        % interfer BS location
        BS_theta = [pi/4,pi/2,3*pi/4,5*pi/4,3*pi/2,7*pi/4];
        BS_location = zeros(2,1,Num_BS);
        for i = 1:Num_BS
            BS_location(:,:,i) = r*[cos(BS_theta(i));sin(BS_theta(i))];  
        end
        
    end







%% first figure, Random User deploy
% draw the circle first
theta = linspace(0, 2*pi, 1000);
% outter circle
x_max = r*cos(theta);
y_max = r*sin(theta);
% second circle
x_second = d(2)*cos(theta);
y_secnod = d(2)*sin(theta);

figure()
plot(x_max,y_max,"-"); % plot the outter circle
hold on;grid on;
plot(center(1),center(2),"b-O","MarkerSize",10);text(center(1),center(2),"BS");
hold on;
plot(x_second,y_secnod,"-","Color",orange);                                % plot the maximum distance that user can deployed
hold on;
plot(User_location(1,:),User_location(2,:),"*","MarkerSize",9);            % user location
text(User_location(1,:),User_location(2,:),["1","2","3","4","5","6",'7','8','9','10']);
xlim([-600 600]);ylim([-600 600]);xlabel("x-axis");ylabel("y-axis");
hold on;
plot(BS_location(1,:),BS_location(2,:),"O",'Color',"r");                   % other BS location
legend('',"Center BS","",'Users',"Interfer BS");title("Random Users Deployment");
text(BS_location(1,:),BS_location(2,:),["BS_1","BS_2","BS_3","BS_4","BS_5","BS_6"]);

