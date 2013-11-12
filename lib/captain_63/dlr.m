function [fit,SE,par,V,COMP,o1,y0] = dlr(y,z,IRW,nvr,alpha,P0,x0,smooth,ALG,p)
% DLR  Dynamic Linear Regression analysis
%
% [fit,fitse,par,parse,comp,e,y0]=dlr(y,z,TVP,nvr,alpha,P0,x0,sm,ALG)
%                                     1 2  3   4    5   6  7  8   9
%
% y: Time series (*)
% z: Regressors (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% nvr: NVR hyper-parameters (0)
% alpha: alpha hyper-parameters for SRW model (1)
% P0: Initial P matrix (1e5)
% x0: Initial state vector (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
% ALG: Smoothing algorithm: P (0) or Q (1-default)
%
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% comp: Linear components
% e: Normalised innovations; use e=e(~isnan(e)) to remove NaNs
% y0: Interpolated data
%
% Example: dlr(y, [ones(size(u)) u], [0 1], 0.001, [0.95 1])
%   regression type model y = c1 + c2(t)*u, with an AR(1) model for
%   c1 (alpha=0.95) and an IRW model for c2 (NVR=0.001 in both cases)
%
% See also DLROPT, FCAST, STAND, DHR, DHROPT, SDP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) and associated m regressors z 
% are specified by the user. Here, z has the same number of rows 
% as y, with a column for each regressor. The function 
% automatically handles missing values in y. In fact, y may be 
% appended with additional NaNs to forecast or backcast beyond 
% the original series. The remaining input arguments are optional.
% 
% TVP is a vector specifying the model associated with each 
% regression parameter, listed in the same order as the columns of 
% z. Choices include a RW/AR(1) model by default (0) or a 
% IRW/SRW model (1). For the case of AR(1) or SRW models, 
% alpha less than unity specifies the additional parameter, while 
% the default value of unity implies a RW or IRW model. For 
% example, a 1st order autoregressive process requires TVP set to 
% zero and 0<alpha<1, where alpha is the AR(1) parameter. 
% Similarly, for a SRW model, TVP is set to unity and 
% 0<alpha<1, where alpha is the smoothing parameter.
% 
% nvr is a vector of NVR hyperparameters for each regressor 
% where, for example, zero (default) implies time invariant 
% parameters. The initial state vector and diagonal of the 
% P-matrix may be specified using x0 and P0, with default values 
% of 0 and 1e5 respectively. FIS may be turned off by changing 
% sm from its default unity to 0. In this case, the model fit and 
% estimated parameters are their filtered values. This speeds up the 
% algorithm and reduces memory usage in cases when smoothing 
% is not required. Finally, either the P (0) or default Q (1) 
% smoothing algorithms are selected using the ALG input 
% argument. Here, the latter is often more robust for RW/IRW 
% models, while SRW models require use of the former. In 
% general, should convergence problems be encountered, changing 
% the algorithm in this manner may help.
% 
% If the lengths of TVP, nvr, alpha, P0 or x0 are less than m, then 
% they are automatically expanded to the correct dimensions by 
% using the final element of the specified input vector. For 
% example, if z has 3 columns but TVP is defined as [1 0], then 
% TVP is automatically expanded to [1 0 0]. Similarly, a scalar P0 
% implies an identity matrix scaled by this value.
% 
% The function returns the model fit (with the same dimensions as 
% y) and parameters par (one column for each regressor), together 
% with the associated standard errors in each case, fitse and parse. 
% 
% It also returns each of the linear components of the model comp, 
% the normalised innovations sequence e and interpolated data y0, 
% where the latter consist of the original series with any missing 
% data replaced by the model. Note that fit is the sum of the 
% columns in comp, while the normalised innovations are padded 
% with initial NaNs to ensure that the vector is the same size as y. 
% If statistical tests on the innovations are required, remove these 
% NaNs with the command e = e(~isnan(e)).

if nargin==0
  disp(' ')
  disp(' DLR  Dynamic Linear Regression analysis')
  disp(' ')
  disp(' [fit,fitse,par,parse,comp,e,y0]=dlr(y,z,TVP,nvr,alpha,P0,x0,sm,ALG)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, z=[]; end
if nargin<3, IRW=[]; end
if nargin<4, nvr=[]; end
if nargin<5, alpha=[]; end
if nargin<6, P0=[]; end
if nargin<7, x0=[]; end
if nargin<8, smooth=[]; end
if nargin<9, ALG=[]; end
if nargin<10, p=[]; end

[fit,SE,par,V,COMP,o1,y0]=dlr0(y,z,IRW,nvr,alpha,P0,x0,smooth,ALG,p);

% end of m-file