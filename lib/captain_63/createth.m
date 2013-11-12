function TH=createth(th_in, P, a, b)
% CREATETH  Modifies THETA matrix
%
% th=createth(th_in,P,a,b)
%
% th_in: Theta matrix (see 'help theta') (*)
% P: Parameter covariance matrix (*)
% b/a: Optional model with the same structure
%      as the original model (trunciated form)
%
% th: Theta matrix (see 'help theta')
%
% See also MCPAR

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Here, th is the THETA matrix with information about the transfer
% function model structure, estimated parameters and their estimated
% accuracy (see help theta), while P is the parameter covariance
% matrix. The 3rd and 4th input arguments are the truncated form
% system numerator b (with assumed unit delay) and denominator a
% (with assumed leading unity) polynomials. If these are omitted
% then the parameter values in th are left unchanged and only the
% covariance matrix component is modified. The function returns a
% modified THETA matrix. It is typically used to evaluate the
% robustness of Proportional-Integral-Plus (PIP) control systems
% using Monte Carlo Simulation.

if nargin==0
  disp(' ')
  disp(' CREATETH  Modifies THETA matrix')
  disp(' ')
  disp(' th=createth(th_in,P,a,b)')
  disp(' ')
  return
end

if nargin<1, th_in=[]; end
if nargin<2, P=[]; end
if nargin<3, a=[]; end
if nargin<4, b=[]; end

TH=createth0(th_in, P, a, b);

% end of m-file