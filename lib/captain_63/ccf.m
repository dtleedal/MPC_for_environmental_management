function out = ccf(y, x, ncoef, s, Title, output)
% CCF  Sample Cross-Correlation Function
%
% tab=ccf(y,x,ncoef,s,title,out)
%
% y: Output time series (*)
% x: Input time series (*)
% ncoef: Number of CCF coefficients (24)
% s: Samples in the season (12)
% title: title for graphical output (automatic)
% out: Tabular and graphical output on (1-default) or off
%
% tab: Tabular output [SCCF, standard errors, Q stats, P-values]
%
% See also ACF, STATIST

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The CCF measures the linear correlation between the selected 
% output and input variables at different lags. It is defined for 
% positive and negative value of the lag. Positive values mean that 
% the output leads the input and vice versa. The function displays a 
% standardised plot of both the time series and the CCF in 
% graphical format (with confidence bands so that a significance 
% test for each coefficient may be completed) and numerical 
% output.
% 
% The time series vectors y (column) and x are the only 
% compulsory inputs to this function.
% 
% The input argument ncoef is an indication of the number of 
% coefficient to calculate. The total number of coefficients is 
% 2*ncoef+1, because it starts at -ncoef and ends at ncoef, 
% including lag 0, where the latter is the contemporaneous linear 
% correlation between both variables.
% 
% Other inputs to this functions are s, the seasonal period that is 
% used in the graphs for a quick location of such values; title is a 
% string variable indicating the required figure title; and out sets 
% the graphical and numerical output on or off, depending on 
% whether the user wishes only to compute the CCF or prefers to 
% visualise the output.
% 
% The output is a table displayed on the MATLAB(r) command 
% window. Here, the columns are the ccf values; their standard 
% deviation; the Ljung-Box Q statistic of cross correlation; and its 
% probability value, i.e. the probability that the Q statistics leaves 
% at the right hand side of its distribution.

if nargin==0
  disp(' ')
  disp(' CCF  Sample Cross-Correlation Function')
  disp(' ')
  disp(' tab=ccf(y,x,ncoef,s,title,out)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, x=[]; end
if nargin<3, ncoef=[]; end
if nargin<4, s=[]; end
if nargin<5, Title=[]; end
if nargin<6, output=[]; end

out=ccf0(y, x, ncoef, s, Title, output);

% end of m-file