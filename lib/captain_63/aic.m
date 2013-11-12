function [arpoly] = aic(x,mod,out)
% AIC  Akaike Information Criterium for Auto-Regressive models
%
% arpoly=aic(y,mod,out)
%
% y: Time series (*)
% mod: Maximum order (scalar) or range of orders (vector 1 x 2)
%      default - min([N/2 32]), where N is series length
% out: Written and graphical output on (1) or off (0-default)
%
% arpoly: AR polynomial of model with minimum AIC
%
% See also MAR, UNIV, UNIVOPT

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is the only compulsory input 
% to this function.
% 
% The second input argument mod defines the range of AR model 
% orders for which to compute the AIC. It may be a scalar, in 
% which case the function will search from order 1 to the specified 
% order; or it may be a vector of dimension 2 indicating the 
% minimum and maximum AR orders to search in. The default 
% value is the minimum of N/2 and 32, where N is the series 
% length. Finally, out specifies whether function should provide 
% tabular and graphical output (1) or not (0 - default).
% 
% The function returns arpoly, i.e. the AR polynomial of the 
% model selected by the AIC.

if nargin==0
  disp(' ')
  disp(' AIC  Akaike Information Criterium for Auto-Regressive models')
  disp(' ')
  disp(' arpoly=aic(y,mod,out)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, mod=[]; end
if nargin<3, out=[]; end

[arpoly]=aic0(x,mod,out);

% end of m-file