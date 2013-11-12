function [J] = IAGP_optimize_SO2(X, A, B, C, D, the_forcing,SO2_emissions,...
    X_hat, Y_hat, optimize_weights, DIdes, DE, g_k_NH, g_k_SH, g_k_SI, offset, phi_bl, PI_bl, use_last_observation,N,i,dt)

%disp(X)

%MPC optimisation test
%AJ - 01-13
% X                     SO2 emissions values to optimize
% A, B, C, D            state space model vectors/matrices
% the_forcing           NH and SH input forcing (already *2)
% X_hat                 state matrix
% Y_hat                 output matrix (NH SH delta temp)
% DI_des                target series
% DE                    integral of error series
% g_k                   present adaptive gain value
% phi_bl                emissions to forcing coefficient
% PI_bl                 NH temp to sea ice extent coefficient
% N                     number of forecast steps
% i                     present time index
% dt                    number of discrete time steps in a year
% DU    the control input sequence
% DIs   the smoothed sea extent sequence
%

%Optimal input pattern
%SO2_emissions = zeros(N,2); % NH and SH SO2 emissions for forecast period
SO2_emissions(i+1:i+N) = X(1:end);
%SO2_emissions(:,1) = X; % put optimization pars in column 1
%SO2_emissions(i+1:i+N) = X(1:end);

% build the forcing to pass to stateSpaceSim

the_forcing(i+1:i+N,1) = the_forcing(i+1:i+N,1) + phi_bl.* SO2_emissions(i+1:i+N)';
%mismatchFactor = 1;

% [X_k, Y_k] = stateSpaceSim(X_k_1,U_k,dt,A,B,C,D)
for j_count = i+1:i+N
    % DIdes is now a timeseries rather than a scalar
    [X_hat_k, Y_hat_k] = stateSpaceSim(...
        X_hat(:,j_count-1),...
        the_forcing(j_count,:)',...
        dt,...
        A,B,C,D);
    X_hat(:,j_count) = X_hat_k; 
    
    Y_hat(:,j_count) = Y_hat_k;
    
    % we need to use the adaptive gain offset (add it on, apply the gain,
    % then subtract it
    adjusted_temp = (Y_hat(1,j_count)+offset).*g_k_NH - offset;
    sea_ice = PI_bl.*adjusted_temp;
    adjusted_sea_ice = (sea_ice + offset).*g_k_SI - offset;
    DE(j_count) = DE(j_count-1) + DIdes(j_count) - adjusted_sea_ice;
    
end
% calculate min sea ice extent and put in x
adjusted_temp = (Y_hat(1,i:i+N)'+offset).*g_k_NH - offset;
    sea_ice = PI_bl.*adjusted_temp;
    x = (sea_ice+offset).*g_k_SI - offset;
% x is the forecast sea ice minimum extent
%disp(x)

%Objective function
% optimize_weights.w1  error between ice and target (100)
% optimize_weights.w2  smoothness of ice (1)
% optimize_weights.w3  integral of error (100)
% optimize_weights.w4  target overshoot changed from zero in 2044 (0)
% optimize_weights.w5  smoothness of emissions (0.05)

%error1 = (DIdes(i+1:i+N-1) - Y_hat(4,i+1:i+N-1)').^2;
% here I am testing whether it is better to include
% the last assimilated sea ice estimate
% if selected we can choose not to use the present observation in the
% definition of the cost this allows pure model inversion
if use_last_observation == 1;
    horizon_error_start = i;
else
    horizon_error_start = i+1;
    x = x(2:end);
end

%Find target overshootovershoot
overshoot_location = find(x>DIdes(horizon_error_start:i+N)); 

error1 = (DIdes(horizon_error_start:i+N) - x).^2;
error2 = var(diff(x));
error3 = (DE(horizon_error_start:i+N).^2);
error4 = (x(overshoot_location)-DIdes((horizon_error_start-1)+overshoot_location)).^2; % brought online in 2044
error5 = sum((diff(SO2_emissions(horizon_error_start+1:i+N,1))).^2);

%TODO error where J is nan if w4 is being used


J = optimize_weights.w1.*sum(error1) + ...
    optimize_weights.w2.*sum(error2) + ...
    optimize_weights.w3.*error3(end) + ...%sum(error3) + ...
    optimize_weights.w4.*sum(error4) + ...
    optimize_weights.w5.*sum(error5);

%disp(J)

%J = sum(w1*(DIdes(i:i+N-1) - DIs(i:i+N-1)).^2) + sum(w2*diff(DU(i:i+N-1)).^2) + w3*sum(DE(i:i+N-1).^2);
% J = sum(w1*(DIdes(i:i+N-1) - DY(i:i+N-1)).^2) + w3*sum(DE(i:i+N-1).^2) ...
%   + sum(w2*diff(DY(i:i+N-1)).^2) + sum(w4*(x(j)-y(j)).^2);

% J = sum(w1*(DIdes(i:i+N-1) - Y_hat(end,i:i+N-1)').^2) + w3*sum(DE(i:i+N-1).^2) ...
%     + sum(w2*diff(Y_hat(end,i:i+N-1)').^2) + sum(w4*(x(j)-y(j)).^2) + w5*var(SO2(i+1:i+N-1));
% w1*((DIdes(i:i+N-1) - Y_hat(end,i:i+N-1)').^2).*[certainty(1) certainty]'
% (DE(i:i+N-1).^2).*[certainty(1) certainty]'
% w2*(diff(Y_hat(end,i:i+N-1)').^2).*certainty'
% w4*((x(j)-y(j)).^2)
% w4*((x(j)-y(j)).^2).*[certainty]'

% J = sum(w1*((DIdes(i:i+N-1) - Y_hat(end,i:i+N-1)').^2).*certainty2') + w3*sum((DE(i:i+N-1).^2).*certainty2') ...
%     + sum(w2*(diff(Y_hat(end,i:i+N-1)').^2).*certainty2(2:end)') + sum(w4*((x(j)-y(j)).^2).*certainty2(j)') + w5*var(SO2(i+1:i+N-1));
% J = sum(w1*((DIdes(i:i+N-1) - Y_hat(end,i:i+N-1)').^2).*certainty2') + w3*sum((DE(i:i+N-1).^2).*certainty2')...
%     + sum(w2*(diff(Y_hat(end,i:i+N-1)').^2).*certainty2(2:end)') + sum(w4*((x(j)-y(j)).^2).*certainty2(j)') + w5*var(SO2(i+1:i+N-1));

% J = sum(w1*(DIdes(i+1:i+N-1) - Y_hat(end,i+1:i+N-1)').^2) + w3*sum((DE(i:i+N-1).^2)) ...
%     + sum(w2*diff(Y_hat(end,i:i+N-1)').^2) + sum(w4*(x(j)-y(j)).^2) + w5*sum(diff(SO2_emissions(i:i+N-1).^2));

%---------------------------------------------------
% original function
% function [J] = MPCoptimisation(X,DU,DF,DIs,DIdes,BSO2,BGHG,A,N);
% %MPC optimisation for IAGP sea ice experiment
% %AJ - 01-13
%
% %Optimal SO2 pattern
% %i = find(X<0);X(i) = 0;
% DU = zeros(N,1);
% DU(1:end-1) = X;
%
% %Desired
% DIdes = ones(N,1)*DIdes;
%
% %Sea ice pattern
% DI = zeros(N,1);DI(1) = DIs;
%
% for i = 1:N-1
%
%    DI(i+1) = BSO2(2)*DU(i) + BGHG(2)*DF(i) - A(2)*DI(i);
%
% end
%
% %Objective function
% w1 = 1;
% w2 = 1e-3;
% J = sum(w1*(DIdes(2:end) - DI(2:end)).^2) + sum(w2*DU(1:end-1).^2);