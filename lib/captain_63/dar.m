function [fit,SE,par,V,COMP,o1,H,y0] = dar(y,p,IRW,nvr,alpha,P0,x0,smooth,ALG,PL,res,Nc)
% DAR  Dynamic Auto-Regression and time frequency analysis
%
% [fit,fitse,par,parse,comp,e,H,y0]=dar(y,na,TVP,nvr,alpha,P0,x0,sm,ALG,PL,res,Nc)
%                                       1 2   3   4    5   6  7  8   9  10 11  12
%
% y: Time series. Standardise first e.g. y=stand(y) (*)
% na: AR order definition (1)
%       Scalar: Full model
%       Vector: Subset AR
%         e.g. na=[1 2 11:14] implies a3,...,a10 = 0
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% nvr: NVR hyper-parameters (0)
% alpha: alpha hyper-parameters for SRW model (1)
% P0: Initial P matrix (1e5)
% x0: Initial state vector (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
% ALG: Smoothing algorithm: P (0) or Q (1-default)
% PL: Plot type (0)
%      -1: No plot and H is empty (saves memory)
%       0: No plot but H is returned
%       1: 3d coloured surface with contours
%       2: 2d contour plot
%       3: Stacked plot
% res: Resolution of H (6)
% Nc: number of contours (PL=2) or stacking distance (PL=3) (10)
%
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% comp: Autoregresive components
% e: Normalised innovations; use e=e(~isnan(e)) to remove NaNs
% H: DAR spectrum e.g. use surf(H)
% y0: Interpolated data
%
% Example: dar(y, [1 12], [1 0], 0.001, [0.95 1])
%   autoregression type model y(k) = a1*y(k-1) + a2*y(k-12) with an SRW 
%   for c1 (alpha=0.95) and an RW model for c2 (NVR=0.001 in both cases)
%
% See also DAROPT, DARSP, FCAST, STAND, DARX, DTFM

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is the only compulsory input 
% to this function. The function automatically handles missing 
% values in y. In fact, y may be appended with additional NaNs to 
% forecast or backcast beyond the original series. It is usually 
% preferable to standardise y before using this function although 
% this is left up to the user, see e.g. 
% help stand. The remaining input arguments are optional. The 
% AR model structure is defined by na, which is a scalar or vector 
% listing the required past output variables used in the model. For 
% example, [1:5, 20] specifies a model based on y(t-1) to y(t-5)
% plus a y(t-20) component (i.e. subset AR).
% 
% TVP is a vector specifying the model associated with each AR 
% model parameter, listed in order of higher powers of the 
% backward shift operator L. Choices include a RW/AR(1) model 
% by default (0) or a IRW/SRW model (1). For the case of AR(1) 
% or SRW models, alpha less than unity specifies the additional 
% parameter, while the default value of unity implies a RW or 
% IRW model. For example, a 1st order autoregressive process 
% requires TVP set to zero and 0<alpha<1, where alpha is the 
% AR(1) parameter. Similarly, for a SRW model, TVP is set to 
% unity and 0<alpha<1, where alpha is the smoothing parameter.
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
% Use PL to specify a time spectra graph. The default zero 
% calculates the spectra and returns the output argument H but 
% does not display the graph, while -1 returns an empty H which 
% saves memory and increases computation time. A value of 1 
% graphs a 3d coloured surface with contours, 2 a 2d contour plot 
% and 3 a stacked plot. The associated input arguments res and Nc 
% control the dimension of H, namely 2^res by length(par), and 
% the number of contours (PL = 2) or stacking distance (PL = 3) 
% respectively. Experiment with res and Nc to obtain the best plot 
% for different computer systems and data sets. The value of Nc is 
% ignored when PL = 0 or PL = 1.
% 
% If the lengths of TVP, nvr, alpha, P0 or x0 are less than the AR 
% model order, then they are automatically expanded to the correct 
% dimensions by using the final element of the specified input 
% vector. For example, if z has 3 columns but TVP is defined as [1 
% 0], then TVP is automatically expanded to [1 0 0]. Similarly, a 
% scalar P0 implies an identity matrix scaled by this value.
% 
% The function returns the model fit (with the same dimensions as 
% y) and parameters par (one column for each regressor), together 
% with the associated standard errors in each case, fitse and parse. 
% It also returns each of the linear components of the model comp, 
% the normalised innovations sequence e and interpolated data y0, 
% where the latter consist of the original series with any missing 
% data replaced by the model. Note that fit is the sum of the 
% columns in comp, while the normalised innovations are padded 
% with initial NaNs to ensure that the vector is the same size as y. 
% If statistical tests on the innovations are required, remove these 
% NaNs with the command e = e(~isnan(e)). Finally, H is the DAR 
% spectrum; for example, use surf(H) to plot the time spectra as a 
% surface.

if nargin==0
  disp(' ')
  disp(' DAR  Dynamic Auto-Regression and time frequency analysis')
  disp(' ')
  disp(' [fit,fitse,par,parse,comp,e,H,y0]=dar(y,na,TVP,nvr,alpha,P0,x0,sm,ALG,PL,res,Nc)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, p=[]; end
if nargin<3, IRW=[]; end
if nargin<4, nvr=[]; end
if nargin<5, alpha=[]; end
if nargin<6, P0=[]; end
if nargin<7, x0=[]; end
if nargin<8, smooth=[]; end
if nargin<9, ALG=[]; end
if nargin<10, PL=[]; end
if nargin<11, res=[]; end
if nargin<12, Nc=[]; end

[fit,SE,par,V,COMP,o1,H,y0]=dar0(y,p,IRW,nvr,alpha,P0,x0,smooth,ALG,PL,res,Nc);

% end of m-file