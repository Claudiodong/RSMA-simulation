% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/06

%% function used to generate the RIS-MU-MISO system channel, with direct and reflect path
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - Nt              (1x1 double) = The number of transmit antenna
% - N_UE            (1x1 double) = The number of user
% - a_d             (1x1 double) = The path loss component for direct path
% - d_r             (1x1 double) = The distance between the RIS and user
% - d_R             (1x1 double) = The distance between the RIS and
%                                   transmitter
% - dis_depe_coeff  (1x1 double) = The distance coefficient coefficient
%                                  which used to model part of path loss
% - a_r             (1x1 double) = The path loss component for reflect path
% - kappa           (1x1 double) = Rician fading factor
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - h_d_k (Nt*N_UE) = The direct channel between transmitter and users
%
% - h_r_k(Nt*N_UE) = The reflect channel between the RIS and users
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function [h_d_k,h_r_k]=Channel_d_r(Nt,N_UE,a_d,d_r,d_R,dis_depe_coeff,a_r,kappa,M)

        %% Direct path
        % MISO system h_direction_d_k(:,1) for user 1
        h_direction_d_k = (randn(Nt,N_UE) + 1i*randn(Nt,N_UE))./sqrt(2);

        % Distance  between user and transmitter
        B_k = pi*rand(1,N_UE) - pi/2; % the Angle B = [-pi/2,pi/2]
        d_k = sqrt( d_r^2 + d_R^2 - 2*d_r*d_R.*cos(B_k));

        % Path loss in linear
        PL = dis_depe_coeff*(d_k.^(-a_d));

        % Direct path 
        % PL(1).*h_direction_d_k(:,1) == h_d_k(:,1) for user 1
        h_d_k = sqrt(PL).*h_direction_d_k;



        %% Reflect path

        % NLOS direction channel following the Rayleight fading 
        % h_direction_k_NLOS(:,1)for user 1
        h_direction_k_NLOS = (randn(M,N_UE)+1i*randn(M,N_UE))./sqrt(2);

        % LOS direction channel      
        h_direction_k_LOS = zeros(M,N_UE); % The direction for th LOS channel
        for i = 1:length(B_k)
            h_direction_k_LOS(:,i) = exp(1i*pi*(0:M-1)*B_k(i));
        end
        % h_direction_k_LOS(:,1) for user 1

        % Rician fading
        % The direction for the Rician channel
        h_direction_r_k = sqrt(kappa/(1+kappa)).*h_direction_k_LOS+sqrt(1/(1+kappa)).*h_direction_k_NLOS;

        % Path loss for the Rician fading
        PL_Rician = dis_depe_coeff*(d_r.^(-a_r));

        % Reflect channel 
        h_r_k = sqrt(PL_Rician).*h_direction_r_k;


%         % plot for user location
%         Tx_location = [0;0]; % x and y
%         RIS_location = [d_R;0]; % RIS is 30 m away the Tx
%         UE_location = RIS_location-[d_r*cos(B_k);d_r*sin(B_k)];
% 
%         figure()
%         plot(Tx_location(1),Tx_location(2),'X');
%         hold on;
%         plot(RIS_location(1),RIS_location(2),'V');
%         hold on;
%         plot(UE_location(1,:),UE_location(2,:),'O')

end