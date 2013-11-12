function sea_ice = IAGP_get_n_step_ice(X, A, B, C, D, the_forcing,...
    X_hat,Y_hat,g_k_NH,g_k_SH, g_k_SI, offset, phi_bl, PI_bl, N,i,dt)

% X                     SO2 emissions values to optimize
% A, B, C, D            state space model vectors/matrices
% the_forcing           NH and SH input forcing (already *2)
% X_hat                 state matrix
% Y_hat                 output matrix (NH SH delta temp)
% g_k                   present adaptive gain value
% phi_bl                emissions to forcing coefficient
% PI_bl                 NH temp to sea ice extent coefficient
% N                     number of forecast steps
% i                     present time index
% dt                    number of discrete time steps in a year
% DU    the control input sequence
% DIs   the smoothed sea extent sequence
%
SO2_emissions = zeros(N,2); % NH and SH SO2 emissions for forecast period
SO2_emissions(:,1) = X; % put optimization pars in column 1
%SO2_emissions(i+1:i+N) = X(1:end);

% build the forcing to pass to stateSpaceSim

the_forcing(i+1:i+N,:) = the_forcing(i+1:i+N,:) + phi_bl.* SO2_emissions;
%mismatchFactor = 1;

% [X_k, Y_k] = stateSpaceSim(X_k_1,U_k,dt,A,B,C,D)
for j_count = i+1:i+N
    
    [X_hat_k, Y_hat_k] = stateSpaceSim(...
        X_hat(:,j_count-1),...
        the_forcing(j_count,:)',...
        dt,...
        A,B,C,D);
    X_hat(:,j_count) = X_hat_k; 
    
    Y_hat(:,j_count) = Y_hat_k;
end
% calculate min sea ice extent and put in sea_ice
adjusted_temp = (Y_hat(1,i+1:i+N)+offset).*g_k_NH - offset;
sea_ice_1 = PI_bl.*adjusted_temp;
sea_ice = ((sea_ice_1+offset).*g_k_SI - offset);
%sea_ice = ((g_k_SI.*(PI_bl.*(((Y_hat(1,i+1:i+N)+10).*g_k_NH)-10))+10)-10);
%sea_ice = g_k_SI.*PI_bl.*(Y_hat(1,i+1:i+N)+g_k_NH); % the forecast sea ice minimum extent