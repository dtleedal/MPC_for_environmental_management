function [NVR, ARp, ARt, sARp, sARt] = univopt(y, par, tmodel, nvr, tar, output, Interv);
% UNIVOPT0  Hyper-parameter and Auto-Regressive (AR) polynomial
%           estimation for Trend + AR type models
%
% [nvr,ARp,ARt,ARpse,ARtse]=univopt(y,p,TVP,nvr0,tar,out,Int)
%
% y: Time series (*)
% p: AR order for Perturbations (vector) (0-none)
% TVP: Integration order for Trend (0-RW/IAR, 1-IRW/DIAR) (1)
% nvr0: Initial NVR hyper-parameter for Trend (1)
% tar: AR order for Trend (vector) (0-none)
% out: Tabular output on (1-default) or off
% Int: Locations of variance intervention points (vector 1/0) (0)
%
% nvr: Estimated NVR hyper-parameter for Trend
% ARp: AR polynomial for Perturbations (no initial unity term)
% ARt: AR polynomial for Trend (no initial unity term)
% ARpse: Standard Errors of AR parameters for Perturbations
% ARtse: Standard Errors of AR parameters for Trend
%
% The variations possible with this function are:
%   Trend only: RW; IRW; IAR; DIAR.
%   Perturbations only: AR.
%   Both: RW+AR; IRW+AR; IAR+AR; DIAR+AR.
%
% Example: [nvr, ARp]=univopt(y, [1:20], 2, 1e-3)
%   IRW+AR. Estimate the optimal NVR for an IRW trend (nvr0=1e-3)
%   together with the polynomial of an AR(20) perturbations model
%
% See also UNIV, AIC, MAR

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% function automatically handles missing values in y. In fact, y 
% may be appended with additional NaNs to forecast or backcast 
% beyond the original series. The remaining input arguments are 
% optional, although most should be normally specified as outputs 
% from the univopt function.
% 
% The first input argument p is a vector indicating the lags that 
% must be used in the auto-regression for the perturbational 
% component about the trend. TVP is a scalar specifying the 
% model associated with the trend. Choices include a RW/AR(1) 
% model by default (0) or a IRW/SRW model (1). nvr0 is the 
% NVR value for the trend. It is normally constrained by the user 
% to avoid identification problems. Although univopt provides 
% default values, in most practical instances these input arguments 
% should be specified by the user, since it is unlikely that the 
% defaults will yield a sensible model for a given data set.
% 
% When a IAR/DIAR trend is required, the trend AR model order 
% is governed by tar which utilizes the same syntax as p above. 
% For IAR trends set TVP to 0, while for DIAR trends set TVP to 
% 1; in both these cases, tar is the vector of indices for the AR 
% trend polynomial. The next input argument, out, specifies 
% tabular display or the results (1) or not (0). Finally, Int allows 
% for sharp (discontinuous) variance intervention for local changes 
% in the trend at the user supplied intervention points. These need 
% to be defined either manually or by some detection method for 
% sharp local changes. Here, Int should take the same dimensions 
% as y, with positive values indicating variance intervention 
% required.
% 
% The function returns nvr, that is the corrected value of the trend 
% NVR based on nvr0 and the final estimate of the innovations 
% variance. ARp and ARt  are the AR polynomials for the 
% perturbations and the trend AR models, when required. ARpse 
% and ARtse are the standard errors of the perturbation and trend 
% AR parameters.

if nargin==0
  disp(' ')
  disp(' UNIVOPT0  Hyper-parameter and Auto-Regressive (AR) polynomial')
  disp('           estimation for Trend + AR type models')
  disp(' ')
  disp(' [nvr,ARp,ARt,ARpse,ARtse]=univopt(y,p,TVP,nvr0,tar,out,Int)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, par=[]; end
if nargin<3, tmodel=[]; end
if nargin<4, nvr=[]; end
if nargin<5, tar=[]; end
if nargin<6, output=[]; end
if nargin<7, Interv=[]; end

[NVR, ARp, ARt, sARp, sARt]=univopt0(y, par, tmodel, nvr, tar, output, Interv);

% end of m-file