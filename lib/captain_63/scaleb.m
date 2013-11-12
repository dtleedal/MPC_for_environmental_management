function bt = scaleb(bt, s)
% SCALEB  Rescale numerator polynomials
%
% bt=scaleb(b,s)
%
% b: Numerator polynomial(s) (see 'help getpar')
% s: Input scaling factors (see 'help prepz')
%
% bt: Rescaled numerator polynomial(s) with the steady
%     state gain adjusted to match the original data
%
% Example: [z, m, s]=prepz([y u], [], [], [], 1);
%   th=riv(z, [2 1 3 0]); [a, b]=getpar(th); bt=scale(b, s)
%   use PREPZ to scale the input (u) to the same numerical
%   range as the output (y), estimate a model using RIV
%   and GETPAR, rescale the numerator polynomial
%
% See also RIV, RIVID, RIVC, RIVCID, GETPAR, PREPZ, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The input argument b consists of the numerator polynomials of 
% a Transfer Function model, while s is the input scaling factors 
% (see 'help prepz'). The output is the rescaled numerator 
% polynomials with the steady state gain adjusted to match the 
% original data.

if nargin==0
  disp(' ')
  disp(' SCALEB  Rescale numerator polynomials')
  disp(' ')
  disp(' bt=scaleb(b,s)')
  disp(' ')
  return
end

if nargin<1, bt=[]; end
if nargin<2, s=[]; end

bt=scaleb0(bt, s);

% end of m-file