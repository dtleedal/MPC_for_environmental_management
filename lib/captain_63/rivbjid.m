function  [th,stats,e,RR,eef,var,Ps,Pn]=rivbjid(Z,nn,sc,flags,nar,gm)
% RIVBJID  Identification of a backward shift MISO transfer function
%
% [th,stats,e,RR,eef,var,Ps,Pn] = rivbjid(Z,nn,sc,flags,nar,gm)
%
% Or for SID object code output (SID Toolbox required)
%   M = rivbjid(Z,nn,sc,flags,nar,gm)
%
% z: I-O data, [Y,U1,...,Unu] where Y and U are column vectors (*)
% nn: Model order search range, 2 by (2+2*nu) matrix (*)
%       [na,nb(1:nu),ntd(1:nu) nc nd] (row 1: search from, row 2: to)
%       na: denominator, nb: numerators, ntd: time delays
%       nc: noise model denominator (AR), nd: numerator (MA)
% sc: Selection criterion (1)
%       1: YIC symmetric; 2: RT2 (I-O only); 3: BIC; 4: EVN
%       5: YIC asymmetric; 6: YIC quadratic; 7: RT2 (full model); 8-AIC
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
% nar: High AR model order used for ARMA modelling with IVARMA method (50)
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
% RR: Table of results for all models searched
% eef: Residuals from noise model (stochastic model residuals)
% var: Variance of residuals: eef if noise model present, else e (SRIV)
% Ps: Covariance matrix for I-O parameter estimates
% Pn: Covariance matrix for ARMA noise model parameter estimates
%
% Example: th=rivbjid([y u], [1 1 1 1 1; 3 3 5 1 1])
%   output y and input u, estimate TF models in the range [1 1 1]
%   to [3 3 5], with first order ARMA(1,1) noise model and listing
%   the results in order of lowest YIC first (default case)
%
% See also RIVBJ, RIVCBJ, RIVCBJID, GETPARBJ, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Differs from original toolbox function RIVID in that ARMA noise
% models are now allowed, i.e. Box-Jenkins (BJ) model. Please see
% supplementary documentation on Upgraded TF Identification and
% Estimation Routines in CAPTAIN (December 2007) for details.

if nargin==0
  disp(' ')
  disp(' RIVBJID  Identification of a backward shift MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,RR,eef,var,Ps,Pn]=rivbjid(Z,nn,sc,flags,nar,gm)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, sc=[]; end
if nargin<4, flags=[]; end
if nargin<5, nar=[]; end
if nargin<6, gm=[]; end

if nargout==1
  th=rivbjid0(Z,nn,sc,flags,nar,gm);
else
  [th,stats,e,RR,eef,var,Ps,Pn]=rivbjid0(Z,nn,sc,flags,nar,gm);
end

% end of m-file