function out = acf(x, ncoef, s, output, param)
% ACF  Sample Autocorrelation and Partial Autocorrelation
%
% tab=acf(y,ncoef,s,out,par)
%
% y: Time series (*)
% ncoef: Number of ACF and PACF coefficients (38)
% s: Samples in the season (12)
% out: Tabular and graphical output on (1-default) or off
% par: Number of parameters when the series is the residuals (0)
%
% tab: [lag, ACF, PACF, standard errors and Q statistic]
%
% See also ARSPEC, CCF, STATIST

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The sample autocorrelation function measures the linear 
% correlation between a time series and several past values. The 
% representation of these coefficients against the lag is called the 
% Autocorrelation Function (ACF). The Partial Autocorrelation 
% Function (PACF) also measures the linear correlation between 
% different lags of a variable, but when all the intermediate lags 
% have been taken into account simultaneously.
%
% The time series vector y (column) is the only compulsory input 
% to this function.
%
% The input argument ncoef is the number of lags to include in the 
% graphical/tabular output. A typcial value is the length of the 
% series divided by four. Input s is the seasonal period of the time 
% series, which adjusts the plot of all the seasonal coefficients in 
% order to make the identification of such lags more 
% straightforward. Input out turns the tabular output on or off and 
% par is the number of parameters in the model, assuming that the 
% series y are effectively the residuals from a previously estimated 
% model. This input is necessary in order to use the correct degrees 
% of freedom of the distribution of the Ljung-Box Q 
% autocorrelation statistic under the null hypothesis.
%
% The output tab is the output table shown in the MATLAB(r) 
% command window when out is set to 1. Here, the columns are 
% the ACF function and its standard error; the Q statistic and its 
% probability value (the probability at the right hand side of the 
% statistic given by its distribution); and the PACF function and 
% standard errors.

if nargin==0
  disp(' ')
  disp(' ACF  Sample Autocorrelation and Partial Autocorrelation')
  disp(' ')
  disp(' tab=acf(y,ncoef,s,out,par)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, ncoef=[]; end
if nargin<3, s=[]; end
if nargin<4, output=[]; end
if nargin<5, param=[]; end

out=acf0(x, ncoef, s, output, param);

% end of m-file