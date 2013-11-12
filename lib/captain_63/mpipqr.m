function [Q, R]=mpipqr(A_mfd, B_mfd, zw, uw, yw)
% MPIPQR  Linear Quadratic weights for multivariable PIP control
%
% [Q,R]=mpipqr(a_mfd,b_mfd,ew,uw,xw)
%
% a_mfd: Matrix Fraction Description (denominator) (*)
% b_mfd: Matrix Fraction Description (numerator) (*)
% ew: vector of integral of error weights (1)
% uw: vector of input weights (1)
% xw: vector of output weights (1)
%
% Q: State weighting matrix
% R: Input weighting matrix
%
% See also MFDFORM, MFD2NMSS, MPIPINIT, PIPLIB, DLQRI

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The first two input arguments are the Matrix Fraction Description,
% for which amfd represents the 'denominator' parameters and bmfd the
% 'numerator' parameters. The final three input arguments are the
% total weights assigned to the integral-of-error state variables ew,
% input state variables uw and output state variables xw. The latter
% three take default values of unity. The function returns the state
% Q and input R weighting matrices.

if nargin==0
  disp(' ')
  disp(' MPIPQR Linear Quadratic weights for multivariable PIP control')
  disp(' ')
  disp(' [Q,R]=mpipqr(a_mfd,b_mfd,ew,uw,xw)')
  disp(' ')
  return
end

if nargin<1, A_mfd=[]; end
if nargin<2, B_mfd=[]; end
if nargin<3, zw=[]; end
if nargin<4, uw=[]; end
if nargin<5, yw=[]; end

[Q, R]=mpipqr0(A_mfd, B_mfd, zw, uw, yw);

% end of m-file