function v=pip(a,b,r,chk);
% PIP  Univariate PIP pole assignment
%
% v=pip(a,b,p,chk)
%
% b/a: System polynomials (trunciated form) (*)
% p: Desired poles (*)
% chk: If positive, does coprimeness test of a and b (0)
%
% v: PIP controller coefficients (see 'help gains')
%
% See also PIPOPT, PIPCOM, PIPCL, PIPLIB, GAINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Truncated form system numerator b (with assumed unit delay) and
% denominator a (with assumed leading unity) polynomials are used
% to determine the PIP control gain vector v for a vector of pole
% positions p. If the length of p is less than the required number
% of poles, then it is automatically expanded by repeating the last
% pole. For example, if p is a scalar, then all the closed-loop
% poles are assigned to this value. If chk is positive, a and b are
% tested for coprimeness (by default, this step is omitted).

if nargin==0
  disp(' ')
  disp(' PIP  Univariate PIP pole assignment')
  disp(' ')
  disp(' v=pip(a,b,p,chk)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, r=[]; end
if nargin<4, chk=[]; end

v=pip0(a,b,r,chk);

% end of m-file