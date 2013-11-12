function [a,b,c,P,delays] = getpar(th);
% GETPAR  Returns transfer function polynomials from theta matrix
%
% [a,b,c,P,d]=getpar(th)
%
% th: theta matrix (see 'help theta')
%
% a: Denominator polynomial
% b: Numerator polynomial
% c: Noise polynomial
% P: Covariance matix
% d: Delays (continuous time only)
%
% Example: th=riv([y u], [2 1 3 0]); [a, b]=getpar(th)
%   estimate a transfer function model using the RIV function
%   then extract the denominator and numerator polynomials
%
% See also RIV, RIVID, RIVC, RIVCID, MAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Returns the denominator a, numerator b and noise c 
% polynomials, together with the noise covariance P matrix and (in 
% the continuous time case) number of delays d. These variables 
% are all extracted from a previously estimated theta-matrix (see help 
% theta). Each row of b represents the numerator polynomial for 
% each input variable.

if nargin==0
  disp(' ')
  disp(' GETPAR  Returns transfer function polynomials from theta matrix')
  disp(' ')
  disp(' [a,b,c,P,d]=getpar(th)')
  disp(' ')
  return
end

if nargin<1, th=[]; end

[a,b,c,P,delays]=getpar0(th);

% end of m-file