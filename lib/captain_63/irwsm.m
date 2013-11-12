function [t,deriv,fitse,trse,derivse,y0,pp,ers,ykk1,er] = irwsm(y,TVP,nvr,Interv,dt,delt)
% IRWSM  Integrated Random Walk smoothing and decimation
%
% [ys,deriv,fitse,trse,derivse,y0,PkN,ers,ykk1,er]=irwsm(y,TVP,nvr,Int,dt,delt)
%
% y: Time series (*)
% TVP: Model type (RW=0, IRW=1, DIRW=2) (1)
% nvr: NVR hyper-parameter (1605*(1/(2*dt))^4)
% Int: Vector of variance intervention points (0)
% dt: decimation rate (1)
%       with dt>1 only the following returned arguments
%       are decimated: ys, deriv, fitse, trse, derivse, y0
% delt: sampling rate (1)
%
% ys: Decimated (or simply smoothed if dt=1) series
% deriv: Derivatives
% fitse: Standard error of model fit (including observation noise)
% trse: Trend (state) standard error
% derivse: standard error of the remaining states (trend derivs.)
% y0: Interpolated data
% PkN: Full smoothed covariance matrix info
% ers: Smoothed error norm. factor
% ykk1: One-step-ahead predictions of the trend
% er: KF error normalisation factor
%
% See also IRWSMOPT, FCAST, STAND, DHR, DHROPT, SDP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% function automatically handles missing values in y. In fact, y 
% may be appended with additional NaNs to forecast or backcast 
% beyond the original series. The remaining input arguments are 
% optional.
% 
% TVP is a vector specifying the model required. Choices include 
% a RW model (0), an IRW model by default (1) or a double 
% integrated random walk model. nvr is a NVR hyperparameter 
% for the model where, for example, zero implies a straight line for 
% the smoothed series. The default valve for nvr is determined by 
% 1605*(1/(2*dt))^4. Int allows for sharp (discontinuous) local 
% changes in the smoothed series at the user supplied intervention 
% points. These need to be defined either manually or by some 
% detection method for sharp local changes. Here, Int should take 
% the same dimensions as y, with positive values indicating 
% variance intervention required. Finally, dt specifies the sampling 
% rate, where the default unity ensures that the sampling rate is the 
% same as for the input series (the function only smoothes the 
% signal), while an integer greater than 1 will decimate the series 
% by an appropriate degree.
% 
% The function returns the smoothed series t (for dt = 1 this 
% variable will have the same dimensions as y) and the derivatives 
% deriv, together with the associated standard errors err. It also 
% returns the filter frequency response filt, the frequency response 
% h and the associated frequency values w. Finally, the 
% interpolated data y0 consist of the original series with any 
% missing data replaced by the model.

if nargin==0
  disp(' ')
  disp(' IRWSM  Integrated Random Walk smoothing and decimation')
  disp(' ')
  disp(' [ys,deriv,fitse,trse,derivse,y0,PkN,ers,ykk1,er]=irwsm(y,TVP,nvr,Int,dt,delt)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, TVP=[]; end
if nargin<3, nvr=[]; end
if nargin<4, Interv=[]; end
if nargin<5, dt=[]; end
if nargin<6, delt=1; end

[t,deriv,fitse,trse,derivse,filt,w,y0,pp,ers,ykk1,er]=irwsm0(y,TVP,nvr,Interv,dt,delt);

% filt: Frequency response of IRW as a low pass filter
% w: Frequency axis for plotting filt values

% end of m-file
