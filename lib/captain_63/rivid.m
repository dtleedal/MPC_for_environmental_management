function [TH,STATS,E,var,Ps,Pc,y0,RR] = rivid(Z,nn,Sc,flags,gui)
% RIVID  Identification of a backward shift MISO transfer function
%
% [th,stats,e,var,Ps,Pc,y0,rr]=rivid(z,nn,sc,flags,g)
%
% z: I-O data, [Y,U1,...,Unu] where Y and Ui are column vectors (*)
% nn: Model order search range, 2 by (2+2*nu) matrix (*)
%       [na,nb(1:nu),nd(1:nu),nc] (row 1: search from, row 2: to)
%       na: denominator, nb: numerators, nd: delays, nc: noise AR
% sc: Selection criterion: 1-YIC symmetric, 2-RT2, 3-AIC, 4-EVN 
%                          5-YIC asymmetric, 6-YIC quadratic (5)
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
%                     0: no stabilisation
%                     1: stabilise filter and instrument generation
%                     2: stabilise filter only (enables estimation
%                         of marginally unstable systems)
%          (7) Yini: Initial conditions (0)
%                     0: original initial y values
%                     1: mean value of initial y values
% g: Postive scalar to activate graphical interface (0)
%
% th: Theta matrix (see 'help theta')
% stats: Statistics [0,YIC,RT2,AIC,S2,o2,EVN,Ybar,RT2AR]
%          RT2: I-O part only, RT2AR: noise model only
%          Total RT2 = RT2+RT2AR*(1-RT2); YIC uses I-O RT2
% e: Model output errors (y=fit+e)
% var: Residual variance
% Ps: Covariance matrix for I-O parameter estimates
% Pc: Covariance matrix for AR parameter estimates
% y0: Interpolated data
% rr: Table of results for all models searched
%
% Example: th=rivid([y u], [1 1 1 0; 3 3 5 0], 2)
%   output y and input u, estimate models in the range [1 1 1 0]
%   to [3 3 5 0] listing the results in order of highest RT2 first
%
% See also RIV, RIVC, RIVCID, GETPAR, PREPZ, SCALEB, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% This function estimates a set of models for a user specified range 
% of orders, selecting the best one according to YIC, RT2, AIC or 
% EVN (see below). The function generates a table of the 20 best 
% models (if this many solutions have converged) according to the 
% selected criterion. It also provides an optional Graphical User 
% Interface in order to make the identification procedure simpler. 
% The function complements (and automatically calls as a 
% subroutine) riv.  For this reason, the input and output arguments 
% or rivid and riv are similar.
% 
% The argument z is a matrix consisting of the output variable in 
% the first column and the input (or inputs) in the remaining k 
% columns, where k is the number of input variables. 
% The model orders are passed by to the function by means of the 
% nn argument. This is a matrix with the format [nafrom nbfrom 
% ndfrom ncfrom; nato nbto ndto ncto], where nafrom and 
% nato are, respectively, the smallest and highest orders (scalars) 
% of the common denominator polynomial; nbfrom and nbto are 
% row vectors (of dimension k) of the minimum and maximum 
% orders for all the numerator polynomials; ndfrom and ndto are 
% vectors (similarly of dimension k) of the minimum and 
% maximum delays for all the inputs in the model; and ncfrom and 
% ncto are the minimum and maximum orders of the AR model 
% for the perturbations (scalar). These are the compulsory input 
% arguments, while the rest are optional and are set to default 
% values when they are not supplied by the user. 
% 
% The number of models to estimate are all the possible 
% combinations of models ranging from the minimum values for 
% each polynomial order to the maximum. For example, if nn is 
% specified as [0 1 2 0; 2 2 3 0], then the model has a single 
% input and the number of models estimated are 12, i.e.
% 
% [0 1 2 0]; [1 1 2 0]; [2 1 2 0]; [0 2 2 0]; 
% [1 2 2 0]; [2 2 2 0]; [0 1 3 0]; [1 1 3 0];
% [2 1 3 0]; [0 2 3 0]; [1 2 3 0]; and [2 2 3 0].
% 
% The input argument flags is a vector of values [Ni Ft Nr Lr Rc 
% Sc] in which: Ni is the number of IV/SRIV iterations (set to 3 by 
% default); Ft specifies the filtering, where 1 indicates prefiltering 
% using a stabilised denominator polynomial, 2 used a stabilised 
% denominator polynomial for both the instruments and the 
% prefiltering (default) and 0 turns the filtering operation off; Nr is 
% the number of RIV iterations (for the co-ordination between the 
% system and noise models); Lr sets the linear regression method 
% to either standard least squares (0) or a SVD/QR algorithm with 
% tolerance equal to Lr (in case collinearity problems are 
% suspected); and Rc switches between the en-block (0 - default) 
% and recursive (1) solutions. If there are any missing (nan) values 
% in the output, the algorithm automatically switches to recursive 
% mode. Finally, Sc specifies the identification criteria (see 
% below), i.e. YIC (1 - default), RT2 (2), AIC (3) or EVN (4). If 
% only some of the values in flags require changing, the rest may 
% be set to -1 for their default values.
% 
% Any positive value for the final input argument g activates a 
% basic Graphical User Interface that makes the selection process 
% easier. It consists of a list of the best models and a plot of the 
% data with the model fit. The user may change the model and the 
% selection criterion by a simple mouse click, and a graphical 
% window shows displays the fit of the highlighted model. When 
% the return key is pressed, the rivid output arguments are 
% returned to the workspace in the normal way.
% 
% The first output argument th provides information about the 
% estimated model in the form of a standard theta matrix (see 'help 
% theta'), from which the estimated polynomials may be recovered 
% using getpar (see 'help getpar').
% 
% The second output argument (stats) is a vector of various 
% criteria useful for the evaluation of the model, i.e. [cond(P), 
% YIC, RT2, AIC, S2, o2, EVN, Ybar, RT2AR]. Here, YIC, 
% RT2 ( ) and AIC are defined by equations (6.46), 6.47) and 
% (4.48). For each of these criteria, the model fit is determined 
% from the input-output part of the model, while the RT2AR term 
% above refers to the RT2 for the noise model. In this regard, the 
% overall fit may be calculated from RT2+RT2AR*(1-RT2). In 
% stats, cond(P) refers to the conditioning of the P covariance 
% matrix, S2 and o2 are the variance of the residuals and output 
% variable respectively, EVN is the log of the average parameter 
% standard errors and, finally, Ybar is the mean of y0 (see below).
%
% The model output errors are stored in e, with variance var; Ps 
% and Pc are the covariance matrix of the system and noise model 
% parameter estimates respectively; y0 recovers the interpolated 
% data, i.e. the missing output observations are replaced by the 
% estimated values; and, finally, rr is the table of the best 20 
% models selected, as it appears in the MATLAB(r) command 
% window.

if nargin==0
  disp(' ')
  disp(' RIVID  Identification of a backward shift MISO transfer function')
  disp(' ')
  disp(' [th,stats,e,var,Ps,Pc,y0,rr]=rivid(z,nn,sc,flags,g)')
  disp(' ')
  return
end

if nargin<1, Z=[]; end
if nargin<2, nn=[]; end
if nargin<3, Sc=[]; end
if nargin<4, flags=[]; end
if nargin<5, gui=[]; end

[TH,STATS,E,var,Ps,Pc,y0,RR] = rivid0(Z,nn,Sc,flags,gui);

% end of m-file