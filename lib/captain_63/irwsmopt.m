function [nvr,opts,separ] = irwsmopt(y,IRW,meth,Interv, opts, ALG);
% IRWSMOPT  Hyper-parameter estimation for IRWSM
%
% [nvr,opts,parse] = irwsmopt(y,TVP,meth,Int,opts,ALG)
%
% y: Time series (*)
% TVP: Model type (RW=0, IRW=1, DIRW=2) (1)
% meth: Estimation method ('ml')
%     'ml': Maximum Likelihood
%     'f#': Sum of squares of the #-step-ahead forecasting errors
% Int: Vector of variance intervention points (0)
% opts: Optimisation options, see help for OPTIMSET (or FOPTIONS for ALG=3)
% ALG: Optimisation algorithm: 0=fminsearch, 1=fminunc, 2=lsqnonlin (not ML) (0)
%        ALG=3 uses older FMINS based optimisation with interrupt button
%
% nvr: Estimated NVR hyper-parameters
% opts: Returned optimisation options. Type 'help foptions' for details
% parse: Standard Error of hyper-parameters (omit to reduce computation time)
%
% See also IRWSM, DHROPT, STAND, FOPTIONS, FMINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% function automatically handles missing values in y. In fact, y 
% may be appended with additional NaNs to forecast or backcast 
% beyond the original series. The remaining input arguments are 
% optional. TVP is a vector specifying the model required.
% 
% Choices include a RW model (0), an IRW model by default (1) 
% or a double integrated random walk model. The 3rd input 
% argument meth selects the estimation method, where the default 
% Maximum Likelihood 'ml' may be replaced by 'f#' to compute 
% the sum of squares of the #-step-ahead forecasting errors. Int 
% allows for sharp (discontinuous) local changes in the parameters 
% at the user supplied intervention points. These need to be 
% defined either manually or by some detection method for sharp 
% local changes. Here, Int should take the same dimensions as y, 
% with positive values indicating variance intervention required. 
% Note that this option does not influence the NVR estimates 
% when meth>0 (the default option).
% 
% Optimisation options may be set using opt (type help foptions 
% for details), while ALG specifies the optimisation algorithm: 
% fmins (0), fminu (1) or leastsq (2). Here, ALG selects between 
% the more efficient gradient search methods of fminu (see help 
% fminu) and the more robust (especially for discontinuous 
% problems) direct search methods of fmins (see help fmins).
% 
% When meth = 'f#' there is the additional option of using leastsq 
% (see help leastsq). Note that if ALG = 2 and meth = 'ml', then 
% an error occurs, since leastsq cannot be used in the Maximum 
% Likelihood case. The Optimisation Toolbox for MATLAB(r) is 
% required to use fminu or leastsq. The function returns the 
% hyperparameter nvr and opts, where the latter provides 
% confirmation of the options utilised, together with the number of 
% function evaluations etc. (type help foptions for details). 
% Finally, parse are the standard errors of the NVR and hyper-
% parameter. Computation time can sometimes be greatly reduced 
% if this 3rd output argument is omitted from the function call.

if nargin==0
  disp(' ')
  disp(' IRWSMOPT  Hyper-parameter estimation for IRWSM')
  disp(' ')
  disp(' [nvr,opts,parse]=irwsmopt(y,TVP,meth,Int,opts,ALG);')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, IRW=1; end
if isempty(IRW), IRW= 1; end
if nargin<3, meth=[]; end
if nargin<4, Interv=[]; end
if nargin<5, opts=[]; end
if nargin<6, ALG=[]; end

if strcmp('ml', lower(meth)), meth= -10; end
if nargout>2
   [nvr,alpha,opts,amp,separ]=dhropt0(y,0,IRW,meth,-2,[],1,[],opts,ALG,[],[],Interv);
else
   [nvr,alpha,opts]=dhropt0(y,0,IRW,meth,-2,[],1,[],opts,ALG,[],[],Interv);
end

% end of m-file