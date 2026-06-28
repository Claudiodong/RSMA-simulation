function [PL_center_BS_user,PL_interfer_BS_UE] = large_scale_fading( ...
    Center_BS_UEs,BS_location,User_location,shadowing_dB,Nt,Num_User,Num_BS)

 %% Large Scale Fading do not change since the UEs are fixed within a Drop time
        
    % center BS to Users
    % Path loss large scale fading in km and dB, Center BS to UE
    PL_dB = 128.1+37.6*log10(Center_BS_UEs./1000)+shadowing_dB*randn(Nt,Num_User);
    % Path loss in linear, Center BS to UEs
    PL_center_BS_user = 10.^((PL_dB)./10);
        
    % other BS to users Path loss 
    interfer_BS_to_UEs = BS_location - User_location; % Every BS disatnce to all users,[BS1-BS6]
    interfer_BS_to_UEs = sqrt(sum(interfer_BS_to_UEs.^2)); % distance
        
    % Path loss in dB
    PL_dB_other_BS_uEs = 128.1+37.6*log10(interfer_BS_to_UEs./1000)+shadowing_dB*randn(Nt,Num_User,Num_BS);
    % Path loss in linear, other BS to UEs
    PL_interfer_BS_UE = 10.^((PL_dB_other_BS_uEs)./10);    
    PL_interfer_BS_UE = squeeze_size(PL_interfer_BS_UE); % squeeze the size if there is diemension of 1
    % PL_interfer_BS_UE = Num_UEs * Num_BS
    % therefore, (:,1) = The Interfere BS 1 to all UEs
    % (1,:) = The All path loss from all interfere BS to UE1.

end