function [fit,fitse,trend,trendse,comp,y0] = univ(y, par, tmodel, nvr, tar, Interv, smooth)
% UNIV  Trend + Auto-Regression (AR) univariate analysis
%
% [fit,fitse,trend,trendse,comp,y0]=univ(y,ARp,TVP,nvr,ARt,Int,sm)
%
% y: Time series (*)
% ARp: AR polynomial for perturbations ([]-none)
% TVP: Integration order for Trend (0-RW/IAR, 1-IRW/DIAR) (1)
% nvr: NVR hyper-parameter for the trend ([]-none)
% ARt: AR polynomial for trend ([]-none)
% Int: Locations of variance intervention points (vector 1/0) (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
%
% fit: Model fit
% fitse: Standard errors of the fit
% trend: Trend
% trendse: Standard errors of the trend
% comp: Total perturbational component
% y0: Interpolated data
%
% The variations possible with this function are:
%   Trend only: RW; IRW; IAR; DIAR
%   Perturbations only: AR
%   Both: RW+AR; IRW+AR; IAR+AR; DIAR+AR
%
% Example: univ(y, ARp, 2, 1e-3)
%   IRW trend (nvr=1e-3) plus an AR model for the
%   perturbations, i.e. the polynomial ARp
%
% See also UNIVOPT, AIC, MAR

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% function automatically handles missing values in y. In fact, y 
% may be appended with additional NaNs to forecast or backcast 
% beyond the original series. The remaining input arguments are 
% optional, although most should be normally specified as outputs 
% from the univopt function.
% 
% ARp and ARt  are the AR polynomials for the perturbations and 
% the trend models (if required), while TVP is a scalar specifying 
% the model associated with the trend. Options for TVP include a 
% RW/AR(1) model by default (0) or a IRW/SRW model (1). nvr 
% is the NVR value for the trend. Int allows for sharp 
% (discontinuous) local changes in the parameters at the user 
% supplied intervention points. These need to be defined either 
% manually or by some detection method for sharp local changes. 
%
% Here, Int should take the same dimensions as y, with positive 
% values indicating variance intervention required. FIS may be 
% turned off by changing sm from its default unity to 0. In this 
% case, the model fit and estimated parameters are their filtered 
% values. This speeds up the algorithm and reduces memory usage 
% in cases when smoothing is not required.
% 
% By selecting appropriate options for the trend and the 
% perturbations, a wide range of overall model structures are 
% available using this function, including trend only models (RW, 
% IRW, IAR, DIAR, etc.), AR models, or various combinations of 
% the two (RW+AR; IRW+AR; IAR+AR; DIAR+AR; etc.).
% 
% The function returns the model fit (with the same dimensions as 
% y), trend and total seasonal component comp, together with the 
% associated standard errors in each case, fitse and trendse. It also 
% returns the interpolated data y0, where the latter consist of the 
% original series with any missing data replaced by the model. 

if nargin==0
  disp(' ')
  disp(' UNIV  Trend + Auto-Regression (AR) univariate analysis')
  disp(' ')
  disp(' [fit,fitse,trend,trendse,comp,y0]=univ(y,ARp,TVP,nvr,ARt,Int,sm)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, par=[]; end
if nargin<3, tmodel=[]; end
if nargin<4, nvr=[]; end
if nargin<5, tar=[]; end
if nargin<6, Interv=[]; end
if nargin<7, smooth=[]; end

[fit,fitse,trend,trendse,comp,y0]=univ0(y, par, tmodel, nvr, tar, Interv, smooth);

% end of m-file