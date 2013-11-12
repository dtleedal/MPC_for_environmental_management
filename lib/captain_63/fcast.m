function y=fcast(x,n)
% FCAST  Prepare data for forecasting and interpolation
%
% x=fcast(y,n)
%
% y: Time series (*)
% n: Select any combination of fore/backcasting or interpolation (*)
%      scalar: forecast n steps ahead
%      row vector: [0 p] Forecast p steps ahead (as above)
%                  [p 0] Backcast p steps
%                  [p q] Interpolate from sample p to sample q
%      matrix with 2 columns: any combination of the above (see example)
%
% x: Time series y but with Not-a-Number (NaN) variables added to
%      indicate to the other modelling functions where forecasting,
%      backcasting and/or interpolation is required
%
% Example: x=fcast(y, [14 19; 22 22; 0 10; 5 0]); fit=dhr(x, 0, 1)
%   Interpolate using DHR over samples 14 to 19 and sample 22 of
%   the original data set, forecast 10 samples beyond the last datum
%   and backcast 5 samples before the start of the time series.
%
% See also STAND, NAN

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column), together with the options n 
% are selected by the user to generate an equivalent time series x 
% but with Not-a-Number (NaN) variables added to indicate to the 
% other CAPTAIN toolbox modelling functions where forecasting, 
% backcasting and/or interpolation is required. For the input 
% argument n, a scalar p or vector [0 p] appends the series with p 
% NaNs, while [p 0] prepends the series with p NaNs for 
% backcasting and [p q] replaces samples p to q with NaNs for 
% interpolation exercises.

if nargin==0
  disp(' ')
  disp(' FCAST  Prepare data for forecasting and interpolation')
  disp(' ')
  disp(' x=fcast(y,n)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, n=[]; end

y=fcast0(x,n);

% end of m-file