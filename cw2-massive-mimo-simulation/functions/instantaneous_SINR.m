

% since the Nr is 1 in this case, the size of the matrix is being reduced
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% - H_center_t     (Q*Nt complex)        = The scheduled user channels from center BS
% - H_interfere_t  (Q*Nt*Num_BS complex) =  The scheduled user channels from interfere BS
% - w_center       (Nt*Q complex)        = ZFBF precoder design for the center BS to users
% - w_interfere    (Q*Nt*Num_BS complex) = random precoder design for the interfere BS
% - Es_per_channel (1x1 double)          = power allocate to the user channel in linear scale
% - N0             (1x1 double)          = noise power in linear scale
%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [SINR]=instantaneous_SINR(H_center_t,H_interfere_t,w_center,w_interfere,Es_per_channel,N0)
             % find the parameters
             [Q,~,Num_BS] = size(H_interfere_t);

             % Compute the interference power from interfere BS
             I_all = zeros(Q,Num_BS);
             for i = 1:Num_BS           
                  I_all(:,i) = diag(Es_per_channel.*abs(H_interfere_t(:,:,i)*w_interfere(:,:,i)).^2);
             end
    
             % Co-scheduled users Power (interference)
             I_all_center_BS = (Es_per_channel.*abs(H_center_t*w_center).^2);
             UE_Power = diag(I_all_center_BS);
             SINR = zeros(1,Q);                      % define the SINR size
             for j = 1:Q                
                % interference from interfere BS to UE
                I_c = sum(I_all(j,:));
                % co-scheduled user interference
                I_scheduled = sum(I_all_center_BS(j,:)) - UE_Power(j);
                % SINR for Scheduled user
                SINR(j) = UE_Power(j)/(I_scheduled+I_c + N0);
                %SINR(j) = UE_Power(j)/(I_c + N0);
             end
end
