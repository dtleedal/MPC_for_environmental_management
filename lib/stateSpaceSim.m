function [X_k, Y_k] = stateSpaceSim(X_k_1,U_k,dt,A,B,C,D)


% all inputs are forcing rather than emissions so
% there is no need to divide any input by number of
% inter-year time steps


X = X_k_1;

for i = 1:1/dt
    X = A*X + B*U_k;% + (Q.*dt)*(noise_sys_k + ar_1*dt*noise_sys_k_1);
    Y = C*X + D*U_k;% + (R.*dt)*(noise_obs_k + ar_1*dt*noise_obs_k_1);
    %             noise_sys_k_1 = noise_sys_k;
    %             noise_obs_k_1 = noise_obs_k;
end
X_k = X;
Y_k = Y;