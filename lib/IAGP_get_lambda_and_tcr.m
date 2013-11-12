function [L, tcr] = IAGP_get_lambda_and_tcr( A,B,C,D,nL,dt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

end_year = 5000;

t = (1:end_year)';
forcing = zeros(size(t));
% 1% per year increase until 2xCO2 (that takes 70 years)
for i = 2:70
    forcing(i) = forcing(i-1) + (3.76 + forcing(i-1))./100;
end

forcing(71:end) = 3.74;
all_X_states = zeros(2*nL,length(t));
all_Y_states = zeros(2,length(t));
for time_index = 2:end_year
    
    
    last_state = all_X_states(:,time_index-1);
    
    
    
    forcings = [forcing(time_index); forcing(time_index)].*2;
    % do Kalman Filter step
    [X_k, Y_k] = stateSpaceSim(last_state,forcings,dt,A,B,C,D);
    
    
    all_X_states(:,time_index) = X_k;
    % place results into matrix
    all_Y_states(:,time_index) = Y_k;
    
    
    
end
L = mean(all_Y_states(:,end));
tcr = mean(all_Y_states(:,70));
%disp(['The LAMBDA is: ' num2str(mean(all_Y_states(:,end)))])


end

