function [th,stats,e,eef,var,Ps,Pn,y0,thd,statsd] =...
  rivcbj(Z,nn,flags,C,nar,gm)

% RIVCBJ  Estimation of a continuous time MISO transfer function
%
% [th,stats,e,eef,var,Ps,Pn,y0,thd,statsd]=rivcbj(Z,nn,flags,C,nar,gm)
%
% Or for SID object code output (SID Toolbox required)
%   M = rivcbj(Z,nn,flags,C,nar,gm)
% 
% z: I-O data, [Y,U1,...,Unu] where Y and Ui are column vectors (*)
% nn: Model order, [na,nb(1:nu),ntd(1:nu) nc nd] (*)
%       na: denominator, nb: numerators, ntd: time delays
%       nc: noise model denominator (AR), nd: numerator (MA)
% flags: Additional parameters vector [Ni,dt,ddt,cf]
%        Missing value or -1 implies use default value in brackets
%          Ni: When Ni>1: Number (maximum) of SRIV iterations (1e-4)
%              (set large to invoke auto-convergence criterion)
%              when Ni is set<1 then it is interpreted as convcrit:
%              convergence criterion to compare with quadratic norm of
%              parameter vector change between iterations (default 1e-4)
%          dt: Sampling interval for continuous time estimation (1)
%          ddt: Sampling for initial discrete time estimation (1)
%          cf: Constant (1) or adaptive (0) pre-filter flag (0)
% C: Prefilter polynomial parameter (1)
%      Polynomial: Polymomial with order 'na'
%      Scalar <=0: Polynomial created automatically from pole value
%                  e.g. zero to utilise multiple integrators, <0 for
%                  stable filter (chosen so that 1/C roughly matches
%                  the bandwidth of the system being estimated)
%      1: Estimate discrete-time filter and convert to continous-time
%         Normally C<=0 is best for fast and C=1 for coarse sampled data
%         (or, if fast sampled, can use C=1 with ddt set at a more coarse
%         subsampling rate for initial discrete-time estimation)
% nar: High AR model order used for ARMA modelling with IVARMA method (50)
% gm: ARMA estimation algorithm (default - 1 when SID present, else 0)
%       0: IVARMA (CAPTAIN)
%       1: PEM (rather quicker that IVARMA but requires SID Toolbox)
%
% th: Theta matrix when called with multiple output arguments ('help theta')
%     If only one output argument: provides model structure of SID type
%    (if SID Toolbox exists) else a structure with the same fields as SID
% stats: Statistics [cond(P),YIC,RT2,BIC,S2,o2,0,Ybar,RT2c,AIC]
%          (RT2: I-O; RT2c: full model)
% e: System model output errors (y=fit+e)
% eef: Residuals from noise model (stochastic model residuals)
% var: Variance of residuals: eef if noise model present, else e (SRIV)
% Ps: Covariance matrix for I-O parameter estimates
% Pn: Covariance matrix for AR parameter estimates
% y0: Interpolated data
% thd: Theta matrix for initial discrete time model
% statsd: Stats vector for discrete time model
%
% Example: [th,stats,e,eef]=rivcbj([y u],[4 3 7 1 1],[-1 0.005],-gamma,[],0);
%   output y and input u, estimate model with a [4 3 7] TF structure and a
%   first order ARMA(1,1) noise model using RIV with default Ni and sampling
%   interval 0.005; starting MMF prefilter with user-defined parameter
%   (-gamma); default AR order; and IVARMA noise model estimation
%
% See also RIVBJ, RIVBJID, RIVCBJID, GETPARBJ, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Differs from original toolbox function RIVC in that ARMA noise
% models are now allowed, i.e. Box-Jenkins (BJ) model. Please see
% supplementary documentation on Upgraded TF Identification and
% Estimation Routines in CAPTAIN (December 2007) for details.

if nargin==0
  disp(' ')
  disp(' RIVCBJ  Estimation of a continuous time MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,eef,var,Ps,Pn,y0,thd,statsd]=rivcbj(Z,nn,flags,C,nar,gm)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, flags=[]; end
if nargin<4, C=[]; end
if nargin<5, nar=[]; end
if nargin<6, gm=[]; end

if nargout==1
  th=rivcbj0(Z,nn,flags,C,nar,gm);
else
  [th,stats,e,eef,var,Ps,Pn,y0,thd,statsd]=rivcbj0(Z,nn,flags,C,nar,gm);
end

% end of m-file