function [a,b,c,d,P,delays]=getparbj(th);
% GETPARBJ  Returns transfer function polynomials from theta matrix
%
% [a,b,c,d,P,T]=getparbj(th)
%
% th: theta matrix (see 'help theta')
%
% a: Denominator polynomial
% b: Numerator polynomial
% c: Noise numerator polynomial
% d: Noise Denominator polynomial
% P: Covariance matix
% T: Delays (continuous time only)
%
% Example: th=rivbj([y u], [2 1 3 0]); [a, b]=getparbj(th)
%   estimate a transfer function model using the RIVBJ function
%   then extract the denominator and numerator polynomials
%
% See also RIVBJ, RIVBJID, RIVCBJ, RIVCBJID, MAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Differs from original toolbox function RIVC in that ARMA noise
% models are now allowed, i.e. Box-Jenkins (BJ) model. Please see
% supplementary documentation on Upgraded TF Identification and
% Estimation Routines in CAPTAIN (December 2007) for details.

if nargin==0
  disp(' ')
  disp(' GETPARBJ  Returns transfer function polynomials from theta matrix')
  disp(' ')
  disp(' [a,b,c,d,P]=getparbj(th)')
  disp(' ')
  return
end

if nargin<1, th=[]; end

[a,b,c,d,P,delays]=getparbj0(th);

% end of m-file