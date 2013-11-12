function [k,p] = dlqri(a,b,q,r,del)
% DLQRI  Iterative linear quadratic regulator design
%
% [k,p]=dlqri(a,b,q,r,del)
%
% a,b: State space form (*)
% q,r: State and input weights (*)
% del: Convergence tolerence (1e-8)
%
% k: Optimal feedback gain matrix
% p: Final P matrix
%
% See also PIPOPT, PIPCOM

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Has similar functionality to dlqr from the MATLAB Control Toolbox.
% Implemented in an iterative form to deal with singular state
% transition matrices. The state transition matrix a, input vector b
% (or matrix in the multivariable case), state weighting matrix q
% and input weighting scalar r (or matrix in the multivariable case)
% are specified by the user. The convergence tolerance del may be
% optionally specified but is usually left at the default value of 1e-8.
% The function returns the converged control gain vector k (or matrix
% in the multivariable case) and P matrix from the discrete time matrix
% Riccati equation.

if nargin==0
  disp(' ')
  disp(' DLQRI  Iterative linear quadratic regulator design')
  disp(' ')
  disp(' [k,p]=dlqri(a,b,q,r,del)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, q=[]; end
if nargin<4, r=[]; end
if nargin<5, del=[]; end

[k,p] = dlqri0(a,b,q,r,del);

% end of m-file