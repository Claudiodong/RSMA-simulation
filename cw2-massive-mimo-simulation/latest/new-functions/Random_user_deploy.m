

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - Input:
% - d_max(1x1 double)= maximum distance from the center Base station that user can be
% placed
%
% - d_min(1x1 double) = minimum distance from the center Base station that user can be
% placed
%
% - Num_User(1x1 double) = the number of user that require in the system
% - center_BS(2x1 double) =  the center BS location
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%
% - User_location (2xNum_User double) =  the cartesian coordinate of users
% relative to center BS
%
%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [User_location,User_distance] = Random_user_deploy(d_max,d_min,Num_User,center_BS)

    % center BS to users
    User_distance = (d_max-d_min)*rand(1,Num_User)+d_min; % random distance in meter
    U_theta = 2*pi*randn(1,Num_User);                     % random  => following the distrubution
    U_x_distance = User_distance.*cos(U_theta)+center_BS(1);        % user x-distance 
    U_y_distance = User_distance.*sin(U_theta)+center_BS(2);        % user y-diatance
    User_location = [U_x_distance;U_y_distance];

end