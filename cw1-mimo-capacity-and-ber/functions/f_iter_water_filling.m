% Claudio Dong MSc CSP Wireless communciation 
% 2024.1.18 Night

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Power Allocation method that using "Water-Filling Algorithm", which is
% iterative algortihm that produce optimal u value (water-level) that
% allocate non-negative power to each channel and obey the power constaint.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  - SNR = signal to noise ratio in decimal (double 1x1)
%
%  - H = Channel (NtxNr, complex) 
%
%  - Pt = transmit power 1x1 double, usually is 1
%
%  - p_tolerance = the power constaint tolerance, to ensure the sum of
% allocated power is 1
%
%  - Nr_MIMO = number of receiver antenna (double 1x1)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Output:
%  - power = the power allocated to each channel based on its strength to
% maximise the data rate
%
%  - Q = the covariance matrix after allocated power on it.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [power,Q] = f_iter_water_filling(SNR,H,Nr_MIMO,Pt,p_tolerance)
        % Perfrom SVD to find the diagonal matrix of signal
        [~,S,V] = svd(H);
        % find the signal power 
        lambda = diag(S.^2);

        % Iterative water filling algorithm (power allocation)
        i = 0;                             % initialise the count
        % loop tp find the constant value miu, that allocate power to all
        % channels with non-negative power.
        while(true)
            i = i+1;                       % +1 for each iteration
            inv = 1./(SNR.*lambda);        % compute the inverse of channel and SNR

            % compute the water level
            u(i) = (1/(Nr_MIMO-i+1))*(1+sum(inv(1:Nr_MIMO-i+1)));

            % allocate power to channels
            power = max(u(i) - inv,0);

            % The 'if' statement is used to determine whether the allocated power is
            % under power constarint, usually it 1. (whether it is converge)
            if ( abs(sum(power)-Pt) < p_tolerance )
                % if the sum of allocated power to abs(channels - Pt) is less
                % than tolerance value, then it break and leave the loop
                break;
            end
        end
        diag_power = zeros(length(V),1);
        % A     ssign the diag_power with allocated power
        diag_power(1:length(power)) = power; 
        % Generate Covariance Matrix
        Q = V*(diag(diag_power))*V';
end