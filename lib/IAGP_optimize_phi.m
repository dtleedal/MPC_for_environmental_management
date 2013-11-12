function [  ] = IAGP_optimize_phi( PI, Faero, volcano, net_TOA, GHG, SO2, DT, SI, SO, L   )
%Called every 10 years to reoptimize SO2 to forcing coefficient and
%NH temp to sea ice coefficient
%   Detailed explanation goes here
%% define optimizer parameters
OPTIONS = optimset('maxfunevals',1e10,'maxiter',1e10);

LB = ones(N,1).*0; % optimization lower bound
UB = ones(N,1).*30; % upper bound
X0 = ones(N,1)*2; % optimization initial parameter values

%% call optimization
X = lsqnonlin('optimize_phi',X0,LB,UB,OPTIONS)

end



function err = optimize_phi(  )

end