function [th,stats,e,eef,vr,Ps,Pn,y0,AH,AHse,PH,Pr,AHn,AHnse]=...
    rivbj(Z,nn,flags,a0,P0,nar,gm)
% RIVBJ  Estimation of a backward shift MISO transfer function
%
% [th,stats,e,eef,var,Ps,Pn,y0,AH,AHse,PH,Pr] = ...
%                                       rivbj(Z,nn,flags,a0,P0,nar,gm)
%                                             1  2   3   4   5  6  7
%
% Or for SID object code output (SID Toolbox required)
%   M = rivbj(z,nn,flags,a0,P0,nar,gm)
%
% z: I-O data, [Y,U1,...,Unu] where Y and U are column vectors (*)
% nn: Model order, [na,nb(1:nu),ntd(1:nu) nc nd] (*)
%       na: denominator, nb: numerators, ntd: time delays
%       nc: noise model denominator (AR), nd: numerator (MA)
% flags: Additional parameters vector [Ni,Ft,Lr,Rc,Stb,Yini]
%        Missing value or -1 implies use default in brackets
%          (1) Ni: When Ni>1: Number (maximum) of SRIV iterations (1e-4)
%                (set large to invoke auto-convergence criterion)
%                when Ni is set<1 then it is interpreted as convcrit:
%                convergence criterion to compare with quadratic norm of
%                parameter vector change between iterations (default 1e-4)
%          (2) Ft: Filtering in IV/SRIV (2)
%                    1: Stabilised A for prefiltering only
%                    2: Stabilised A for instruments and prefilter
%                    0: Filtering turned off
%          (3) Lr: Linear regression method (0)
%                    0: Standard
%                    tol>0: SVD/QR robust algorithm
%          (4) Rc: Block (0-default) or recursive algorithm (1-slower)
%          (5) Stb: Stabilisation of filter and model polynomial (1)
%                     0: No stabilisation
%                     1: Stabilise filter and instrument generation
%                     2: Stabilise filter only (enables estimation
%                        of marginally unstable systems)
%          (6) Yini: Initial conditions (0)
%                     0: Original initial y values
%                     1: Mean value of initial y values
% a0: Initial parameter estimates for recursive algorithm (0)
% P0: Initial covariance diagonal for recursive algorithm (1000)
% nar: High AR model order used for ARMA modelling with IVARMA method (20)
% gm: ARMA estimation algorithm (default - 1 when SID present, else 0)
%       0: IVARMA (CAPTAIN)
%       1: PEM (rather quicker that IVARMA but requires SID Toolbox)
%
% th: Theta matrix when called with multiple output arguments ('help theta')
%     If only one output argument: provides model structure of SID type
%     (if SID Toolbox exists) else a structure with the same fields as SID
% stats: Statistics [cond(P),YIC,RT2,BIC,S2,o2,EVN,Ybar,RT2c,AIC,YICa,YICt]
%        (RT2: I-O; RT2c: full model; YICa: asymmetric P; YICt: quadratic)
% e: System model output errors (y=fit+e)
% eef: Residuals from noise model (stochastic model residuals)
% var: Variance of residuals: eef if noise model present, else e (SRIV)
% Ps: Covariance matrix for I-O parameter estimates
% Pn: Covariance matrix for AR parameter estimates
% y0: Interpolated data
% AH: Recursive estimation parameters history (with recursive solution only)
% AHse: Recursive estimates of standard deviations of AH (ditto)
% PH: Recursive estimates covariance matrix evolution
%       e.g. AH(:,24) is the parameters vector estimate at sample 24
%       with the standard deviations AHse(:,24) and the covariance matrix
%       given by PH(:,23*size(AH,1)+(1:size(AH,1)))
% Pr: Asymmetric RIV covariance matrix estimate
%
% Example: [th, stats]=rivbj([y u], [2 1 3 1 1], 1e-5, [], [], [], 1)
% output y and input u, estimate model with a [2 1 3] TF structure,
% using convergence criterion norm 1e-5 (or 20 IV iterations if less)
% with ARMA(1,1) model estimated by PEM gradient algorithm if available
%
% See also RIVBJID, RIVCBJ, RIVCBJID, GETPARBJ, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Differs from original toolbox function RIV in that ARMA noise
% models are now allowed, i.e. Box-Jenkins (BJ) model. Please see
% supplementary documentation on Upgraded TF Identification and
% Estimation Routines in CAPTAIN (December 2007) for details.

if nargin==0
  disp(' ')
  disp(' RIVBJ  Estimation of a backward shift MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,eef,vr,Ps,Pn,y0,AH,AHse,PH,Pr,AHn,AHnse]=rivbj(Z,nn,flags,a0,P0,nar,gm)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, flags=[]; end
if nargin<4, a0=[]; end
if nargin<5, P0=[]; end
if nargin<6, nar=[]; end
if nargin<7, gm=[]; end

if nargout==1
  th=rivbj0(Z,nn,flags,a0,P0,nar,gm);
else
  [th,stats,e,eef,vr,Ps,Pn,y0,AH,AHse,PH,Pr,AHn,AHnse]=rivbj0(Z,nn,flags,a0,P0,nar,gm);
end

% end of m-file
