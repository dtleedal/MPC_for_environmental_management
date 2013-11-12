function out = statist(x,output)
% STATIST  Sample descriptive statistics
%
% tab=statist(y,out)
%
% y: Time series vector or matrix of columns (*)
% out: Table output on (1-default) or off
%
% tab: Tabular output returned as a variable
%
% See also ACF, CCF

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% This function returns various statistics for each column vector of 
% variables in the matrix y. The list of statistics includes: the 
% number of samples; number of missing observations; minimum; 
% first quartile; median; third quartile; maximum; mean; geometric 
% mean; range; interquartile range; standard deviation; variance; 
% mean/standard deviation; skewness; kurtosis; and the Jarque-
% Bera test.
% 
% The only compulsory input to this function is y. The second 
% input arguments (out) turns the on screen display of these 
% statistics on (1 - default) or off (0).
% 
% The output argument tab is simply the table shown in the 
% command window, returned as a matrix.

if nargin==0
  disp(' ')
  disp(' STATIST  Sample descriptive statistics')
  disp(' ')
  disp(' tab=statist(y,out)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, output=[]; end

out=statist0(x, output);

% end of m-file