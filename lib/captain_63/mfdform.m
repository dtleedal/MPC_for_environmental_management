function [A,B] = mfdform(a,b)
% MFDFORM  Matrix Fraction Description form
%
% [a_mfd,b_mfd]=mfdform(a,b)
% 
% a: Denominator matrix polynomial, formed row by row without the
%      intial unity and padded with zeros to uniform length (*)
% b: Numerator matrix polynomial, formed row by row, padded with
%      zeros and with unit delay assumed (*)
%
% a_mfd: Matrix Fraction Description (denominator)
% b_mfd: Matrix Fraction Description (numerator)
%
% See also MFD2NMSS, MPIPQR, MPIPINIT, PIPLIB, DLQRI

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Returns the MFD form, for which amfd represents the ‘denominator’
% parameters and bmfd the ‘numerator’ parameters. The denominator
% matrix polynomial a is formed row by row without the initial
% unity and is padded with zeros to uniform length. The numerator
% matrix polynomial b is also formed row by row, padded with zeros
% and with unit delay assumed. The format required for a and b is
% illustrated by the on-line multivariable Proportional-Integral-
% Plus (PIP) control demonstration PIPDEMO3.

if nargin==0
  disp(' ')
  disp(' MFDFORM  Matrix Fraction Description form')
  disp(' ')
  disp(' [a_mfd,b_mfd]=mfdform(a,b)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end

[A,B]=mfdform0(a,b);

% end of m-file