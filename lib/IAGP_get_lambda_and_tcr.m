function [L, tcr] = IAGP_get_lambda_and_tcr( A,B,C,D,nL,dt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

end_year = 5000;

t = (1:end_year)';
forcing = zeros(size(t));

pre_industrial_CO2 = 278;
CO2_conc = pre_industrial_CO2;

% make the forcing series for a 1% per year CO2 increase to pre-ind doubling
for i = 1:100
    forcing(i) = 5.35*log(CO2_conc/278);
    CO2_conc = CO2_conc + 0.01*CO2_conc;
    if (CO2_conc > 2*pre_industrial_CO2); CO2_conc = 2*pre_industrial_CO2; end
end

forcing(101:end) = forcing(100);


all_X_states = zeros(2*nL,length(t));
all_Y_states = zeros(2,length(t));
for time_index = 2:end_year
    
    
    last_state = all_X_states(:,time_index-1);
    
    
    
    forcings = [forcing(time_index); forcing(time_index)];
    % do Kalman Filter step
    [X_k, Y_k] = stateSpaceSim(last_state,forcings,dt,A,B,C,D);
    
    
    all_X_states(:,time_index) = X_k;
    % place results into matrix
    all_Y_states(:,time_index) = Y_k;
    
    
    
end
L = 3.7083/mean(all_Y_states(:,end));
tcr = mean(all_Y_states(:,70));
%disp(['The LAMBDA is: ' num2str(mean(all_Y_states(:,end)))])


end

