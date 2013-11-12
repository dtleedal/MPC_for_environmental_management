function [tfs,fit,SE,par,V,COMP,o1,y0] = darx(y,u,nn,IRW,nvr,alpha,P0,x0,smooth,ALG)
% DARX  Dynamic AutoRegressive multi-eXogenous variables analysis
%
% [tfs,fit,fitse,par,parse,comp,e,y0]=darx(y,u,nn,TVP,nvr,alpha,P0,x0,sm,ALG)
%                                          1 2 3   4   5    6   7  8  9  10
%
% y: Time series (*)
% u: Input (*)
% nn: Model structure [na,nb(1:nu),nd(1:nu)] ([1 1 0])
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% nvr: NVR hyper-parameters (0)
% alpha: alpha hyper-parameters for SRW model (1)
% P0: Initial P matrix (1e5)
% x0: Initial state vector (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
% ALG: Smoothing algorithm: P (0) or Q (1-default)
%
% tfs: Transfer function (simulation) output
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% comp: Model components for y(k-1), ..., u(k-1), ...
% e: Normalised innovations; use e=e(~isnan(e)) to remove NaNs
% y0: Interpolated data
%
% Example: darx(y, [u1 u2], [1 1 1 2 3], 0, 0.001)
%   transfer function type model y(k) = a(k)*y(k-1) + b1(k)*u1(k-2) + b2(k)*u2(k-3)
%   with RW models for both parameters (NVR=0.001)
%
% See also DARXOPT, FCAST, STAND, DAR, DTFM

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for DARX: Paul McKenna

% The time series vector y (column) and input matrix u are 
% specified by the user. Each column of u represents an input 
% signal. The function automatically handles missing values in y. 
% In fact, y may be appended with additional NaNs to forecast or 
% backcast beyond the original series, as long as appropriate 
% values for u are also specified. The remaining input arguments 
% are optional. The DARX model structure is defined by nn, 
% which takes the form [n,m,d] where, in transfer function terms, 
% n and m are the number of denominator and numerator 
% parameters respectively, while d is the number of samples time 
% delay. A first order model with unity time delay and one 
% numerator parameter [1, 1, 1] is utilised by default.
% 
% TVP is a vector specifying the model associated with each 
% DARX model parameter, listed in order of each denominator 
% parameter and then the numerator parameters for each input, 
% i.e. a(1,t)...a(n,t), b(0,t)...b(m,t) for the single input,
% single output example. Choices include a RW/AR(1) model by 
% default (0) or a IRW/SRW model (1). For the case of AR(1) or 
% SRW models, alpha less than unity specifies the additional 
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
% If the lengths of TVP, nvr, alpha, P0 or x0 are less than the 
% total number of parameters, then they are automatically 
% expanded to the correct dimensions by using the final element of 
% the specified input vector. For example, if the DARX model has 
% 3 parameters but TVP is defined as [1 0], then TVP is 
% automatically expanded to [1 0 0]. Similarly, a scalar P0 implies 
% an identity matrix scaled by this value.
% 
% The function returns the simulation response tfs (with the same 
% dimensions as y), regression fit and parameters par (one column 
% for each parameter), together with the associated standard errors 
% in the latter two cases, fitse and parse. Here, tfs is based on 
% feeding the input signal through the model (the output signal is 
% not used, except to establish the initial conditions), while fit 
% represents the 1-step ahead predictions and is equivalent to the 
% fit returned by dlr.
% 
% The function also returns each of the linear components of the 
% model comp, i.e. the components associated with each input and 
% output and their past values, the normalised innovations 
% sequence e and interpolated data y0, where the latter consist of 
% the original series with any missing data replaced by the model. 
% Note that fit is the sum of the columns in comp, while the 
% normalised innovations are padded with initial NaNs to ensure 
% that the vector is the same size as y. If statistical tests on the 
% innovations are required, remove these NaNs with the command 
% e = e(~isnan(e)).

if nargin==0
  disp(' ')
  disp(' DARX  Dynamic AutoRegressive multi-eXogenous variables analysis')
  disp(' ')
  disp(' [tfs,fit,fitse,par,parse,comp,e,y0]=darx(y,u,nn,TVP,nvr,alpha,P0,x0,sm,ALG)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, u=[]; end
if nargin<3, nn=[]; end
if nargin<4, IRW=[]; end
if nargin<5, nvr=[]; end
if nargin<6, alpha=[]; end
if nargin<7, P0=[]; end
if nargin<8, x0=[]; end
if nargin<9, smooth=[]; end
if nargin<10, ALG=[]; end

[tfs,fit,SE,par,V,COMP,o1,y0]=darx0(y,u,nn,IRW,nvr,alpha,P0,x0,smooth,ALG);

% end of m-file