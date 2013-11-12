function [th,stats,e,thd,statsd] = rivc(Z,nn,flags,C,convcrit);
% RIVC  Estimation of a continuous time MISO transfer function
%
% [th,stats,e,thd,statsd]=rivc(z,nn,flags,c,conv)
%
% z: I-O data, [Y,U1,...,Unu] where Y and Ui are column vectors (*)
% nn: Model order, [na,nb(1:nu),nd(1:nu)] (*)
%       na: denominator, nb: numerators, nd: delays
% flags: Additional parameters vector [Ni,dt,ddt,cf]
%        Missing value or -1 implies use default value in brackets
%          Ni: Number (maximum) of SRIV iterations (20)
%          (set large to invoke auto-convergence criterion)
%          dt: Sampling interval for continuous time estimation (1)
%          ddt: Sampling for initial discrete time identification (1)
%          cf: Constant (1) or adaptive (0) pre-filter flag (0)
% C: Prefilter polynomial parameter (1)
%      (a) Polymomial with order 'na'
%      (b) Scalar <=0: Polynomial created automatically from pole value
%       e.g. zero to utilise multiple integrators, <0 for stable filter 
%      (chosen so that 1/C roughly matches the bandwidth of the system
%       being estimated)
%      (c) 1: estimate discrete-time filter and convert to continous-time
%       Normally (b) is best for fast and (c) for coarse sampled data
% conv: convergence criterion to compare with quadratic norm of
%       parameter vector change between iterations (default 1e-4)
%
% th: Theta matrix (see 'help theta')
% stats: Statistics [cond(P),YIC,RT2,-,S2,o2,-,Ybar,-]
% e: Model output errors (y=fit+e)
% thd: Theta matrix for initial discrete time model
% statsd: Stats vector for discrete time model
%
% Example: th=rivc([y u], [2 1 3], [2 -1 -1 1])
%   output y and input u, estimate model with a [2 1 3] structure
%   using 2 SRIV iterations and a constant estimated prefilter
%
% See also RIV, RIVID, RIVCID, GETPAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The argument z is a matrix consisting of the output variable in 
% the first column and the input (or inputs) in the remaining k 
% columns, where k is the number of input variables.
% The model orders are passed by to the function by means of the 
% nn argument. This is a vector with the format [na nb nd], where 
% na is the order (a scalar) of the common denominator 
% polynomial; nb is a vector (of dimension k) of orders for all the 
% numerator polynomials; and nd is a vector (similarly of 
% dimension k) of delays for all the inputs in the model. These are 
% the compulsory input arguments, while the rest are optional and 
% are set to default values when they are not supplied by the user.
% 
% The input argument flags is a vector of values [Ni dt ddt cf] in 
% which: Ni is the number of SRIV iterations (set to 3 by default); 
% dt is the sampling interval for continuous time estimation (1); 
% ddt is the sampling interval for initial discrete time 
% identification; and cf is specifies either constant (1 - default) or 
% adaptive (0) pre-filtering.
% 
% The third input argument C selects the prefiltering polynomial. 
% In the default case (1),  a discrete time filter is estimated and 
% converted to continuous time. Setting C to zero or a negative 
% scalar uses a stable first order filter created automatically from 
% this pole value. This is equivalent to Young's method of multiple 
% filters or the Poisson-Moment-Function (PMF) method. Finally, 
% to specify the filter directly, set C to the required polynomial (a 
% vector of length na).
% 
% The first output argument th provides information about the 
% estimated model in the form of a standard theta matrix (see 'help 
% theta'), from which the estimated polynomials may be recovered 
% using getpar (see 'help getpar').
% 
% The second output argument (stats) is a vector of various 
% criteria useful for the evaluation of the model, i.e. [cond(P), 
% YIC, RT2, 0, S2, o2, 0, Ybar]. Here, YIC and RT2 are 
% defined by equations (6.46) and (4.48), while cond(P) refers to 
% the conditioning of the P covariance matrix, S2 and o2 are the 
% variance of the residuals and output variable respectively and 
% Ybar is the mean of the output. Note that some elements of 
% stats are zero, for compatibility with the equivalent output 
% argument of the function riv.
% 
% The model output errors are stored in e, while thd is the theta 
% matrix (see 'help theta') for the initial discrete time model and 
% statsd is the 'stats' vector for the discrete time model (see 'help 
% riv' for a definition of this latter stats vector).

if nargin==0
  disp(' ')
  disp(' RIVC  Estimation of a continuous time MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,thd,statsd]=rivc(z,nn,flags,c)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, flags=[]; end
if nargin<4, C=[]; end
if nargin<5, convcrit=[]; end

[th,stats,e,thd,statsd] = rivc0(Z,nn,flags,C,convcrit);

% end of m-file