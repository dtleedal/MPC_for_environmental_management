function [per, f] = period(x, opt, vec)
% PERIOD  Estimate periodogram
%
% [per,f]=period(y,out,vec)
%
% y: Time series (*)
% out: Table and graphical output ([1 0])
%        out(1): output on (1) or off
%        out(2): log scale on (1) or off
% v: Vector of periods where vertical lines should be plotted
%
% per: Estimated periodogram
% f: Frequency axis at which the periodogram is estimated
%
% See also ARSPEC

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% optional argument out controls the tabular and graphical output, 
% where out(1) turns the output on (1 by default) or off (0) and 
% out(2) specifies a log (1 by default) or normal scale (0).
% With regards to the output arguments, per represents the 
% estimated periodogram, while f is the frequency axis at which 
% the periodogram is estimated.

if nargin==0
  disp(' ')
  disp(' PERIOD  Estimate periodogram')
  disp(' ')
  disp(' [per,f]=period(y,out,vec)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, opt=[]; end
if nargin<3, vec=[]; end

[per, f]= period0(x, opt, vec);

% end of m-file