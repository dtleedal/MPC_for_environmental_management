function [X_hat_k, Y_hat_k,P_k,Y_var] = kalmanFilterSim(X_k_1,Y_k,U_k,dt,A,B,C,D,Q,R,P,assimilate)
% make the zero order hold

u = U_k.*[1;1;dt];
X = X_k_1;


for i = 1:((1/dt))
    
    X = A*X + B*u;
    
end
P = A*P*A' + Q; % a priori covariance
if (assimilate == 1)
    
    K = (P*C')/(C*P*C' + R);
    
    % calculate Kalman gain
    
    X = X + K*(Y_k - (C*X+D*u));
    % a posteriori covariance
    P = (eye(length(X)) - K*C)*P;
    
    
    % data assimilation or not
    
    %     X = A*X + B*u(:,1/dt);
    %     P = A*P*A' + (dt.*Q); % a priori covariance
end

X_hat_k = X;
Y_hat_k = C*X + D*u;
P_k = P;
Y_var = R + C*P*C';