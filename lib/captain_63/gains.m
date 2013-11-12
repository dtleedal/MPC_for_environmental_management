function [f, g, k, r]=gains(a, b, v)
% GAINS  Univariate PIP control polynomials
%
% [f,g,k,r]=gains(a,b,v)
%
% b/a: System polynomials (trunciated form) (*)
% v: PIP controller coefficients (*)
%
% f: Output polynomial
% g: Input polynomial
% k: Integral gain
% r: Command polynomial (if present)
%
% See also PIP, PIPOPT, PIPCOM, PIPCL, PIPLIB

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Truncated form system numerator b (with assumed unit delay) and
% denominator a (with assumed leading unity) polynomials are used
% to determine the order of the PIP control polynomials. Here, v
% is the control gain vector (e.g. returned by pip or pipopt),
% while f, g and k are the output feedback polynomial, denominator
% polynomial for the input filter and integral gain respectively.
% If v is obtained using pipcom, then r is the optional command
% input polynomial.

if nargin==0
  disp(' ')
  disp(' GAINS  Univariate PIP control polynomials')
  disp(' ')
  disp(' [f,g,k,r]=gains(a,b,v)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, v=[]; end

[f,g,k,r]=gains0(a,b,v);

% end of m-file