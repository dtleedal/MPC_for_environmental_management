function [v,F,g,d,h,Q,R]=pipopt(a,b,erw,uw,xw)
% PIPOPT  Univariate PIP Linear Quadratic optimal control
%
% [v,F,g,d,h,Q,r]=pipopt(a,b,ew,uw,xw)
%
% b/a: System polynomials (trunciated form) (*)
% ew: Control error weight (*)
% uw: Control input weights (*)
% xw: State weights (*)
%
% v: PIP controller coefficients (see 'help gains')
% F: State transition matrix
% g: Input vector
% d: Command vector
% h: Observation vector
% Q: State weighting matrix
% r: Input weighting
%
% See also PIP, PIPCOM, PIPCL, PIPLIB, GAINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Truncated form system numerator b (with assumed unit delay) and
% denominator a (with assumed leading unity) polynomials are used
% to determine the control gain vector v based on PIP-LQ design,
% where ew, uw and xw are the total weights for the integral-of-
% error, input and output state variables. If required, the
% function also returns the non-minimal state space form state
% transition matrix F, input vector g, command input vector d
% and observation vector h, together with the state weighting
% matrix Q and scalar input weight r.

if nargin==0
  disp(' ')
  disp(' PIPOPT  Univariate PIP Linear Quadratic optimal control')
  disp(' ')
  disp(' [v,F,g,d,h,Q,R]=pipopt(a,b,erw,uw,xw)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, erw=[]; end
if nargin<4, uw=[]; end
if nargin<5, xw=[]; end

[v,F,g,d,h,Q,R]=pipopt0(a,b,erw,uw,xw);

% end of m-file