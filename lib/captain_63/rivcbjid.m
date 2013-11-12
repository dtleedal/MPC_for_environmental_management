function  [th,stats,e,RR,eef,var,Ps,Pn]=rivcbjid(Z,nn,sc,flags,C,nar,gm)

% RIVCBJID  Identification of a continuous time MISO transfer function
%
% [th,stats,e,RR,eef,var,Ps,Pn]=rivcbjid(z,nn,sc,flags,C,nar,gm)
%
% Or for SID object code output (SID Toolbox required)
%   M = rivcbjid(Z,nn,flags,C,nar,gm)
%
% z: I-O data, [Y,U1,...,Unu] where Y and U are column vectors (*)
% nn: Model order search range, 2 by (2+2*nu) matrix (*)
%       [na,nb(1:nu),ntd(1:nu) nc nd] (row 1: search from, row 2: to)
%       na: denominator, nb: numerators, ntd: time delays
%       nc: noise model denominator (AR), nd: numerator (MA)
% sc: Selection criterion (1-YIC, 2-RT2, 3-BIC, 4-AIC) (1)
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
% RR: Table of results for all models searched
% eef: Residuals from noise model (stochastic model residuals)
% var: Variance of residuals: eef if noise model present, else e (SRIV)
% Ps: Covariance matrix for I-O parameter estimates
% Pn: Covariance matrix for ARMA noise model parameter estimates
%
% Example: th=rivcbjid([y u], [1 1 1 1 1; 3 3 5 1 1])
%   output y and input u, estimate TF models in the range [1 1 1]
%   to [3 3 5], with first order ARMA(1,1) noise model and listing
%   the results in order of lowest YIC first (default case)
%
% See also RIVBJ, RIVBJID, RIVCBJ, GETPARBJ, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Differs from original toolbox function RIVCID in that ARMA noise
% models are now allowed, i.e. Box-Jenkins (BJ) model. Please see
% supplementary documentation on Upgraded TF Identification and
% Estimation Routines in CAPTAIN (December 2007) for details.

if nargin==0
  disp(' ')
  disp(' RIVCBJID  Identification of a continuous time MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,RR,eef,var,Ps,Pn]= rivcbjid(Z,nn,sc,flags,C,nar,gm)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, sc=[]; end
if nargin<4, flags=[]; end
if nargin<5, C=[]; end
if nargin<6, nar=[]; end
if nargin<7, gm=[]; end

if nargout==1
  th=rivcbjid0(Z,nn,sc,flags,C,nar,gm);
else
  [th,stats,e,RR,eef,var,Ps,Pn]= rivcbjid0(Z,nn,sc,flags,C,nar,gm);
end

% end of m-file
