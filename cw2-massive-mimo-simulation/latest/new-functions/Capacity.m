function C = Capacity(center_H,H_interfere,w_center,w_interfere,S,N0)

         [Q,~,Num_BS]=size(H_interfere);
         % Random precoder for the interfere BSs        
         I_Interfere_all = zeros(Q,Num_BS);

         for i_bs = 1:Num_BS
             %I_Interfere_all(:,i_bs) = diag(S.*abs(H_interfere(:,:,i_bs)*w_interfere(:,:,i_bs)).^2);
             % do not take the digonal value, which consider the sum of all
             % products. for example, all the product of H_1 with all
             % precoder for scheduled user in cell. therefore, it have
             % interference from Q, which is sum all of them.
             L = H_interfere(:,:,i_bs)*w_interfere(:,:,i_bs);
             I_Interfere_all(:,i_bs) = sum(S.*abs(L).^2,2);
         end

         % Compute the SINR for the User
         UE_all_power = S.*abs(center_H*w_center).^2;
         UE_power = diag(UE_all_power);
         SINR = zeros(1,Q);
         C = zeros(1,Q);
         for i_e = 1:Q
             % inter-user interference (Co-scheduled), should be close to 0
             % by ZFBF
             I_c = sum(UE_all_power(i_e,:)) - UE_power(i_e);

             % interference from all 6 interfere BS to UE
             I_interfere = sum(I_Interfere_all(i_e,:));

             % compute the SINR
             SINR(i_e) = UE_power(i_e)/(N0+I_c+I_interfere);

             % Compute the Rate
             C(i_e) = log2(1+SINR(i_e));
         end
end