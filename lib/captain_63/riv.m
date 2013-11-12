function [TH,STATS,E,vr,Ps,P,y0,AH,AHse,PH,Pr] = riv(Z,nn,flags,a0,P0)
% RIV  Estimation of a backward shift MISO transfer function
%
% [th,stats,e,var,Ps,Pc,y0,AH,AHse,PH,Pr]=riv(z,nn,flags,a0,P0)
%
% z: I-O data, [Y,U1,...,Unu] where Y and U are column vectors (*)
% nn: Model order, [na,nb(1:nu),nd(1:nu),nc] (*)
%      (na: denominator, nb: numerators, nd: delays, nc: noise AR)
% flags: Additional parameters vector [Ni,Ft,Nr,Lr,Rc,Stb,Yini]
%        Missing value or -1 implies use default in brackets
%        Ni and Nr are mutually exclusive: there is automatic
%        switching between SRIV (Ni>0, Nr=0; nc=0) and full RIV
%        (Ni=0; Nr>0; nc>0) when necessary
%          (1) Ni: Number of basic IV/SRIV iterations (3 - SRIV}
%          (2) Ft: Filtering in IV/SRIV (2)
%                    1: Stabilised A for prefiltering only
%                    2: Stabilised A for instruments and prefilter
%                    0: Filtering turned off
%          (3) Nr: Number of RIV iterations including
%                    AR noise model (0 if nc=0, else 3)
%          (4) Lr: Linear regression method (0)
%                    0: Standard
%                    tol>0: SVD/QR robust algorithm
%          (5) Rc: Block (0-default) or recursive algorithm (1-slower)
%          (6) Stb: Stabilisation of filter and model polynomial (1)
%                     0: No stabilisation
%                     1: Stabilise filter and instrument generation
%                     2: Stabilise filter only (enables estimation
%                        of marginally unstable systems)
%          (7) Yini: Initial conditions (0)
%                     0: Original initial y values
%                     1: Mean value of initial y values
% a0: Initial parameter estimates for recursive algorithm (0)
% P0: Initial covariance diagonal for recursive algorithm (1000)
%
% th: Theta matrix (see 'help theta')
% stats: Statistics [cond(P),YIC,RT2,AIC,S2,o2,EVN,Ybar,RT2AR,YICa,YICt]
%          RT2: I-O part only, RT2AR: noise model only
%          Total RT2 = RT2+RT2AR*(1-RT2); YIC uses I-O RT2
%          YICa - YIC using asymmetric matrix
%          YICt - YIC 'total' using full quadratic form, not just trace
% e: Model output errors (y=fit+e)
% var: Residual variance
% Ps: Covariance matrix for I-O parameter estimates
% Pc: Covariance matrix for AR parameter estimates
% y0: Interpolated data
% AH: Recursive estimation parameters history (recursive solution only)
% AHse: Recursive estimates of standard deviations of AH 
% PH: Recursive estimates covariance matrix evolution
%       e.g. AH(:,24) is the parameters vector estimate at sample 24
%       with the standard deviations AHse(:,24) and the covariance matrix
%       given by PH(:,23*size(AH,1)+(1:size(AH,1)))
% Pr: Asymmetric RIV covariance matrix estimate
%
% Example: th=riv([y u], [2 1 3 0], [-1 0])
%   output y and input u, estimate a model with a [2 1 3 0]
%   structure, using 3 IV iterations and no filtering
%
% See also RIVID, RIVC, RIVCID, GETPAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The argument z is a matrix consisting of the output variable in 
% the first column and the input (or inputs) in the remaining k 
% columns, where k is the number of input variables.
% The model orders are passed by to the function by means of the 
% nn argument. This is a vector with the format [na nb nd nc], 
% where na is the order (a scalar) of the common denominator 
% polynomial; nb is a vector (of dimension k) of orders for all the 
% numerator polynomials; nd is a vector (similarly of dimension k) 
% of delays for all the inputs in the model; and nc is the order of 
% the AR model for the perturbations (scalar). These are the 
% compulsory input arguments, while the rest are optional and are 
% set to default values when they are not supplied by the user.
% 
% The input argument flags is a vector of values [Ni Ft Nr Lr Rc] 
% in which: Ni is the number of IV/SRIV iterations (set to 3 by 
% default); Ft specifies the filtering, where 1 indicates prefiltering 
% using a stabilised denominator polynomial, 2 used a stabilised 
% denominator polynomial for both the instruments and the 
% prefiltering (default) and 0 turns the filtering operation off; Nr is 
% the number of RIV iterations (for the co-ordination between the 
% system and noise models); Lr sets the linear regression method 
% to either standard least squares (0) or a SVD/QR algorithm with 
% tolerance equal to Lr (in case collinearity problems are 
% suspected); and Rc switches between the en-block (0 - default) 
% and recursive (1) solutions.
% 
% If there are any missing (nan) values in the output, the algorithm 
% automatically switches to recursive mode. If only some of the 
% values in flags require changing, the rest may be set to -1 
% (default values). The final input arguments (a0 and P0) are the 
% initial conditions for the mean and covariance of the parameters. 
% These arguments are ignored in en block mode.
% 
% The first output argument th provides information about the 
% estimated model in the form of a standard theta matrix (see 'help 
% theta'), from which the estimated polynomials may be recovered 
% using getpar (see 'help getpar').
% 
% The second output argument (stats) is a vector of various 
% criteria useful for the evaluation of the model, i.e. [cond(P), 
% YIC, RT2, AIC, S2, o2, EVN, Ybar, RT2AR]. Here, YIC, 
% RT2 and AIC are defined by equations (6.46), 6.47) and 
% (4.48). For each of these criteria, the model fit is determined 
% from the input-output part of the model, while the RT2AR term 
% above refers to the RT2 for the noise model. In this regard, the 
% overall fit may be calculated from RT2+RT2AR*(1-RT2). In 
% stats, cond(P) refers to the conditioning of the P covariance 
% matrix, S2 and o2 are the variance of the residuals and output 
% variable respectively, EVN is the log of the average parameter 
% standard errors and, finally, Ybar is the mean of y0 (see below).
%
% The model output errors are stored in e, with variance var. Ps 
% and Pc are the covariance matrix of the system and noise model 
% parameter estimates respectively. Finally, y0 recovers the 
% interpolated data, i.e. the missing output observations are 
% replaced by the estimated values.

if nargin==0
  disp(' ')
  disp(' RIV  Estimation of a backward shift MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,var,Ps,Pc,y0,AH,AHse,PH,Pr]=riv(z,nn,flags,a0,P0)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, flags=[]; end
if nargin<4, a0=[]; end
if nargin<5, P0=[]; end

[TH,STATS,E,vr,Ps,P,y0,AH,AHse,PH,Pr]=riv0(Z,nn,flags,a0,P0);

% end of m-file