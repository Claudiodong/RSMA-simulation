% Claudio Dong Imperial College Londom
% Wireless Communications and Optimisations 
% Coursework 3
% 2024/03/06

%% function used to generate the RIS-MU-MISO system channel, between the transmitter and RIS
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - Nt              (1x1 double) = The number of transmit antenna
% - N_UE            (1x1 double) = The number of user
% - M               (1x1 double) = The number of element in the RIS
% - d_R             (1x1 double) = The distance between the RIS and
%                                   transmitter
% - dis_depe_coeff  (1x1 double) = The distance coefficient coefficient
%                                  which used to model part of path loss
% - a_r             (1x1 double) = The path loss component for reflect path
% - kappa           (1x1 double) = Rician fading factor
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% - G               (M*Nt complex) = The channel between the transmitter
% and RIS, each antenna to M elements, G(1,:) for first antenna to all
% elements in the RIS.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function [G]=G_channel(Nt,N_UE,M,dis_depe_coeff,d_R,a_r,kappa)
       % channel between the transmitter and RIS
       PL_G = dis_depe_coeff*d_R^(-a_r); % Path loss between RIS and transmitter

       % NLOS for G direction
       G_direction_NLOS = (randn(Nt,M)+1i*randn(Nt,M))./sqrt(2);

       % LOS for G direction
       w_r = (pi/2)*randn(1,Nt);
       w_t = (pi/2)*randn(1,Nt); 
       % should have size of N*1 for each user
       a_rr = zeros(N_UE,N_UE);
       % should have size of M*1 for each user
       a_t = zeros(M,N_UE);

       % The G_LOS should have size N*M, where G_LOS(1,:) for first antenna
       for i = 1:Nt
           % arr should have size (N_UE*1)
           a_rr(:,i) = exp(1i*pi*(0:Nt-1)*w_r(i));
           a_t(:,i) = exp(1i*pi*(0:M-1)*w_t(i));     
       end
       % LOS direction channel
       G_LOS = a_rr(:,1)*a_t(:,1)';
       %G_LOS = a_rr*a_t';
       % G direction
       G_direction = sqrt(kappa/(1+kappa))*G_LOS + sqrt(1/(1+kappa))*G_direction_NLOS;
       % Final G channel from transmitter to RIS
       G = sqrt(PL_G).*G_direction;
       G = G.'; % transpose
end
