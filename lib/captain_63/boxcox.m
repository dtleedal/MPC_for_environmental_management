function xt = boxcox(x, trans, bin, out)
% BOXCOX  Optimal Box-Cox transformation for homoskedasticity
%
% xt=boxcox(x,trans,bin,out)
%
% x: Time series (*)
% trans: Vector of tansformations ([-1:0.5:1])
% bin: Number of divisions for rank-mean and SE-mean graphs
%      (if vector, points at which divisions should appear)
% out: Output off (0) or on (1 - default)
%
% xt: Transformed variable
%
% See also HISTON, CUSUM

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% This function produces several diagnostics in order to look for 
% time varying mean and variance in a time series. It computes 
% several such Box-Cox transformations and assesses each of 
% them in terms of several criteria, including the standard error of 
% the transformed variable (it should be a minimum); the mean 
% over the standard error (maximum); and the likelihood function 
% (maximum).
% 
% In addition, it provides two graphical plots, useful to check for 
% stationarity in the mean and variance simultaneously. In order to 
% build these graphs, the time series is divided into several bins 
% and the range, standard deviation and mean of each bin is 
% estimated. The plots consist of representing the mean against the 
% standard deviation of each bin (standard error-mean plot) and 
% the range against the mean (range-mean plot).
% 
% The time series vector y (column) is the only compulsory input 
% to this function.
% 
% The input trans is a vector of real values that selects the 
% particular transformations in order to compute the criteria with; 
% bin is the number of bins to divide the variable or, if a vector, 
% the limits of each bin; finally, out sets the output on or off.
% This function returns the variable transformed according to the 
% likelihood criteria.

if nargin==0
  disp(' ')
  disp(' BOXCOX  Optimal Box-Cox transformation for homoskedasticity')
  disp(' ')
  disp(' xt=boxcox(x,trans,bin,out)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, trans=[]; end
if nargin<3, bin=[]; end
if nargin<4, out=[]; end

xt = boxcox0(x, trans, bin, out);

% end of m-file