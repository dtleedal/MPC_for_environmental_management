function [fit,fitse,Tend,Tendse,DHR,e,AMP,PHASE,sDHR,sDHRse,y0,dhrse,Xkk,Pkk,ers,ykk1] = dhr(y,P,IRWharm,NVR,alpha,P0,x0,smooth,ALG,Interv,intD)
% DHR  Dynamic Harmonic Regression analysis
%
% [fit,fitse,tr,trse,comp,e,amp,phs,ts,tsse,y0,dhrse,Xhat,Phat,ers,ykk1] ...
%     = dhr0(y,P,IRWharm,NVR,alpha,P0,x0,smooth,ALG,Int,IntD)
%            1 2    3     4    5   6  7    8     9  10   11
%
% y: Time series (*)
% P: Periodic components; set P(1)=0 to include a trend (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW, 2-Trigonommetric) (0)   
%        (for LLT use RW and IRW trends simultaneously)
% nvr: NVR hyper-parameters (0)
% alpha: alpha parameters; set alpha=1 for RW or IRW model (1)
% P0: Initial P matrix (1e6)
% x0: Initial state vector (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
% ALG: Smoothing algorithm: P (0) or Q (1-default)
% Int: Vector of variance intervention points (zeros(length(y),1))
% IntD: Variance intervention matrix diagonal (1e2 for trend level)
%
% fit: Model fit
% fitse: Standard errors of the fit
% tr: Trend (when trend is specified as IRW, second column is trend slope)
% trse: Standard errors of the trend (as above: slope standard errors)
% comp: Harmonic components
% e: Normalised innovations; use e=e(~isnan(e)) to remove NaNs
% amp: Amplitude of harmonic components
% phs: Phase of harmonic components
% ts: Total seasonal component
% tsse: Standard errors of the total seasonal component
% y0: Interpolated data
% dhrse: Standard errors of all the components (not grouped)
% Xhat: Parameters X(k/k) when just KF is run, X(k/N) when both KF/FIS
% Phat: P-matrix P(k/k) when just KF is run, P(k/N) when both KF/FIS
% ers:  Normalised error as returned by kalmsmo: 1+H(k)'*P(k/N)*H(k)
% ykk1: KF one step ahead predictions
%
% Example: dhr(y, [0 12./(1:6)], [1 0], [0.001 0.01], [0.95 1])
%   SRW trend model (NVR=0.001, alpha=0.95) together with 6 periodic 
%   components (12 and harmonics) each modelled with a RW (NVR=0.01)
%
% See also DHROPT, FCAST, STAND, DLR, SDP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) and associated m periodic 
% components P are specified by the user. Set the first element of 
% P to include a trend. For example, [0 12] implies a trend and a 
% seasonal component for monthly data. The function 
% automatically handles missing values in y. In fact, y may be 
% appended with additional NaNs to forecast or backcast beyond 
% the original series. The remaining input arguments are optional.
% 
% TVP is a vector specifying the model associated with each 
% regression parameter, listed in the same order as the elements of 
% P. Choices include a RW/AR(1) model by default (0) or a 
% IRW/SRW model (1). For the case of AR(1) or SRW models, 
% alpha<1 specifies the additional parameter, while the 
% default value of unity implies a RW or IRW model. For 
% example, a 1st order autoregressive process requires TVP set to 
% zero and 0<alpha<1, where alpha is the AR(1) parameter. 
% Similarly, for a SRW model, TVP is set to unity and 
% 0<alpha<1, where alpha is the smoothing parameter. Finally, a 
% LLT model is obtained by using RW and IRW trends 
% simultaneously, i.e. with P set to [0 0].
% 
% nvr is a vector of NVR hyperparameters for each regressor 
% where, for example, zero (default) implies time invariant 
% parameters. The initial state vector and diagonal of the 
% P-matrix may be specified using x0 and P0, with default values 
% of 0 and 1e6 respectively. FIS may be turned off by changing 
% sm from its default unity to 0. In this case, the model fit and 
% estimated parameters are their filtered values. This speeds up the 
% algorithm and reduces memory usage in cases when smoothing 
% is not required. Also, either the P (0) or default Q (1) smoothing 
% algorithms are selected using the ALG input argument. Here, 
% the latter is often more robust for RW/IRW models, while SRW 
% models require use of the former. In general, should 
% convergence problems be encountered, changing the algorithm 
% in this manner may help. Int allows for sharp (discontinuous)
% local changes in the parameters at the user supplied intervention
% points. These need to be defined either manually or by some
% detection method for sharp local changes. Here, Int should take
% the same dimensions as y, with positive values indicating variance
% intervention required. Finally, InD gives the diagonal of the
% variance intervention matrix.
% 
% If the lengths of TVP, nvr, alpha, P0 or x0 are less than the 
% length of P, then they are automatically expanded to the correct 
% dimensions by using the final element of the specified input 
% vector. For example, if P has 3 elements but TVP is defined as 
% [1 0], then TVP is automatically expanded to [1 0 0]. Similarly, 
% a scalar P0 implies an identity matrix scaled by this value.
% 
% The function returns the model fit (with the same dimensions as 
% y), trend tr and total seasonal component ts (i.e. the sum of all 
% the seasonal components, which are returned individually as a 
% matrix comp), together with the associated standard errors in 
% each case, fitse and trse and tsse. It also returns the normalised 
% innovations sequence e, amplitude amp and phase phs of the 
% harmonic components and, finally, the interpolated data y0, 
% where the latter consist of the original series with any missing 
% data replaced by the model. Note that the normalised 
% innovations are padded with initial NaNs to ensure that the 
% vector is the same size as y. If statistical tests on these are 
% required, remove the NaNs with the command e = e(~isnan(e)).

if nargin==0
  disp(' ')
  disp(' DHR  Dynamic Harmonic Regression analysis')
  disp(' ')
  disp(' [fit,fitse,tr,trse,comp,e,amp,phs,ts,tsse,y0,dhrse]=dhr(y,P,TVP,nvr,alpha,P0,x0,sm,ALG,Int,IntD)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, P=[]; end
if nargin<3, IRWharm=[]; end
if nargin<4, NVR=[]; end
if nargin<5, alpha=[]; end
if nargin<6, P0=[]; end
if nargin<7, x0=[]; end
if nargin<8, smooth=[]; end
if nargin<9, ALG=[]; end
if nargin<10, Interv=[]; end
if nargin<11, intD=[]; end

[fit,fitse,Tend,Tendse,DHR,e,AMP,PHASE,sDHR,sDHRse,y0,dhrse,Xkk,Pkk,ers,ykk1]=...
   dhr0(y,P,IRWharm,NVR,alpha,P0,x0,smooth,ALG,Interv,intD);

% end of m-file
