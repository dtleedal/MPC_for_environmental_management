function [TH,STATS,E,RR] = rivcid(Z,nn,flags,C);
% RIVCID  Identification of a continuous time MISO transfer function
%
% [th,stats,e,rr]=rivcid(z,nn,flags,c)
%
% z: I-O data,  [Y,U1,...,Unu] where Y and U are column vectors (*)
% nn: Model order search range, 2 by (2+2*nu) matrix (*)
%       [na,nb(1:nu),nd(1:nu)] (row 1: search from, row 2: to)
%       na: denominator, nb: numerators, nd: delays
% flags: Additional parameters vector [Ni,dt,ddt,cf]
%        Missing value or -1 implies use default value in brackets
%          Ni: Number of SRIV iterations (3)
%          dt: Sampling interval for continuous time estimation (1)
%          ddt: Sampling for initial discrete time identification (1)
%          cf: Constant (1) or adaptive (0) pre-filter flag (0)
%          Sc: Selection criterion (1-YIC, 2-RT2, 3-AIC, 4-EVN) (1)
% C: Prefilter polynomial parameter (1)
%      (a) Polymomial with order 'na'
%      (b) Scalar <=0: Polynomial created automatically from pole value
%       e.g. zero to utilise multiple integrators, <0 for stable filter 
%      (chosen so that 1/C roughly matches the bandwidth of the system
%       being estimated)
%      (c) 1: estimate discrete-time filter and convert to continous-time
%       Normally (b) is best for fast and (c) for coarse sampled data
%
% th: Theta matrix (see 'help theta')
% stats: Statistics [cond(P),YIC,RT2,-,S2,o2,-,Ybar,-]
% e: Model output errors (y=fit+e)
% rr: Table of all the models searched
%
% Example: th=rivcid([y u], [1 1 1; 3 3 5], [-1 -1 -1 -1 2])
%   output y and input u, estimate models in the range [1 1 1]
%   to [3 3 5] listing the results in order of highest RT2 first
%
% See also RIV, RIVID, RIVC, GETPAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% This function estimates a set of models for a user specified range 
% of orders, selecting the best one according to YIC, RT2, AIC or 
% EVN (see below). The function generates a table of the 20 best 
% models (if this many solutions have converged) according to the 
% selected criterion. It also provides an optional Graphical User 
% Interface in order to make the identification procedure simpler. 
% The function complements (and automatically calls as a 
% subroutine) rivc.  For this reason, the input and output 
% arguments or rivcid and rivc are similar.
% 
% The argument z is a matrix consisting of the output variable in 
% the first column and the input (or inputs) in the remaining k 
% columns, where k is the number of input variables.
% This is a matrix with the format [nafrom nbfrom ndfrom; nato 
% nbto ndto], where nafrom and nato are, respectively, the 
% smallest and highest orders (scalars) of the common 
% denominator polynomial; nbfrom and nbto are row vectors (of 
% dimension k) of the minimum and maximum orders for all the 
% numerator polynomials; ndfrom and ndto are vectors (similarly 
% of dimension k) of the minimum and maximum delays for all the 
% inputs in the model. These are the compulsory input arguments, 
% while the rest are optional and are set to default values when 
% they are not supplied by the user.
% 
% The input argument flags is a vector of values [Ni dt ddt cf Sc] 
% in which: Ni is the number of SRIV iterations (set to 3 by 
% default); dt is the sampling interval for continuous time 
% estimation (1); ddt is the sampling interval for initial discrete 
% time identification; and cf is specifies either constant (1 - 
% default) or adaptive (0) pre-filtering. Finally, Sc specifies the 
% identification criteria (see below), i.e. YIC (1 - default),   (2) , 
% AIC (3) or  EVN (4). If only some of the values in flags require 
% changing, the rest may be set to -1 for their default values.
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
% the conditioning of the   covariance matrix, S2 and o2 are the 
% variance of the residuals and output variable respectively and 
% Ybar is the mean of the output. Note that some elements of 
% stats are zero, for compatibility with the equivalent output 
% argument of the function rivid.
% 
% The model output errors are stored in e, while thd is the theta 
% matrix (see 'help theta') for the initial discrete time model and 
% statsd is the 'stats' vector for the discrete time model (see 'help 
% riv' for a definition of this latter stats vector). Finally, rr is the 
% table of the best 20 models selected, as it appears in the 
% MATLAB(r) command window.

if nargin==0
  disp(' ')
  disp(' RIVCID  Identification of a continuous time MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,rr]=rivcid(z,nn,flags,c)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, flags=[]; end
if nargin<4, C=[]; end

[TH,STATS,E,RR] = rivcid0(Z,nn,flags,C);

% end of m-file