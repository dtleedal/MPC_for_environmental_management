function [NVR,ALPHA,opts,separ] = dlropt(y,z,IRW,method,nvr,alpha,nvr0,alpha0,opts,ALG,output,P0,p);
% DLROPT  Hyper-parameter estimation for DLR
%
% [nvr,alpha,opts,parse]=dlropt(y,z,TVP,meth,nvrc,alphac,nvr0,alpha0,opts,ALG,tab,P0)
%                               1 2  3   4    5     6     7     8     9   10  11  12
%
% y: Time series (*)
% z: Regressors (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% meth: Estimation method ('ml')
%     'ml': Maximum Likelihood
%     'f#': Sum of squares of the #-step-ahead forecasting errors
% nvrc: Constraints for each NVR (-2)
%       -2: Free estimation
%       -1: Constrained estimation (all parameters with nvrc=-1 are equal)
%      >=0: NVR constrained to this value (it is not estimated)
% alphac: Constraints for each alpha (-2, -1, or >=0 as for nvrc) (1)
% nvr0: Initial NVR hyper-parameters (0.0001)
% alpha0: Initial alpha hyper-parameters (1)
% opts: Optimisation options, see help for OPTIMSET (or FOPTIONS for ALG=3)
% ALG: Optimisation algorithm: 0=fminsearch, 1=fminunc, 2=lsqnonlin (not ML) (0)
%        ALG=3 uses older FMINS based optimisation with interrupt button
% tab: Display: 0=none, 1=tabulate results, 2=update window (2)
% P0: Initial P matrix (1e5)
%
% nvr: Estimated NVR hyper-parameters
% alpha: Estimated alpha hyper-parameters
% opts: Returned optimisation options. Type 'help foptions' for details
% parse: Standard Error of hyper-parameters (omit to reduce computation time)
%
% Example: dlropt(y, [ones(size(u)) u], 0, [], [0 -2])
%   regression type model y = c1 + c2(t)*u, with an RW model for
%   both parameters and c1 assumed constant (NVR fixed at zero)
%
% See also DLR, FCAST, STAND, FOPTIONS, FMINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) and associated m regressors z 
% are specified by the user. Here, z has the same number of rows 
% as y, with a column for each regressor. The function 
% automatically handles missing values in y. In fact, y may be 
% appended with additional NaNs to forecast or backcast beyond 
% the original series. The remaining input arguments are optional.
% TVP is a vector specifying the model associated with each 
% regression parameter, listed in the same order as the columns of 
% z. Choices include a RW/AR(1) model by default (0) or a 
% IRW/SRW model (1). meth is the estimation method, where the 
% default Maximum Likelihood 'ml' may be replaced by 'f#' to 
% compute the sum of squares of the #-step-ahead forecasting 
% errors.
% 
% nvrc defines the constraints for each NVR, where -2 implies 
% free estimation, -1 constrained estimation (all parameters with 
% nvrc = -1 are equal) and >=0 implies the associated NVR is 
% constrained to this value (it is not estimated). alphac defines 
% similar constraints for each alpha parameter (-2, -1, or >=0 as for 
% nvrc). Initial NVR and alpha hyper-parameters may be specified 
% using nvr0 and  alpha0 respectively. For example, to optimise 
% for a RW or IRW model, ensure alphac and alpha0 are unity 
% (the default). To optimise component 'i' for a SRW, set the i'th 
% element of TVP to unity and alphac(i) to -2, -1, or fixed at 
% 0<alphac<1. This normally produces an improved fit to the 
% spectrum, but computation time is longer.
% 
% Optimisation options may be set using opt (type help foptions 
% for details), while ALG specifies the optimisation algorithm: 
% fmins (0), fminu (1) or leastsq (2). Here, ALG selects between 
% the more efficient gradient search methods of fminu (see help 
% fminu) and the more robust (especially for discontinuous 
% problems) direct search methods of fmins (see help fmins). 
% When meth = 'f#' there is the additional option of using leastsq 
% (see help leastsq). Note that if ALG = 2 and meth = 'ml', then 
% an error occurs, since leastsq cannot be used in the Maximum 
% Likelihood case. The Optimisation Toolbox for MATLAB(r) is 
% required to use fminu or leastsq.
% 
% Finally, tab defines the display options and P0 specifies the 
% initial diagonal of the P matrix. Here, if tab = 1 or 2, then the 
% final results are displayed in tabular form. Additionally, if tab = 
% 2, a window appears during optimisation showing the latest 
% value of the Likelihood Function or the Sum-of-Squares for the 
% #-step-ahead forecasting errors. When ALG = 0 (fmins) and tab 
% = 2, a 'stop button' will appear below the update window: click 
% to terminate the optimisation and return the current estimates.
% If the lengths of TVP, nvrc, alphac, nvr0, alpha0 or P0 or x0 
% are less than m, then they are automatically expanded to the 
% correct dimensions by using the final element of the specified 
% input vector. For example, if z has 3 columns but TVP is 
% defined as [1 0], then TVP is automatically expanded to [1 0 0]. 
% Similarly, a scalar P0 implies an identity matrix scaled by this 
% value.
% 
% The function returns vectors of NVR and alpha hyperparameters, 
% nvr and alpha respectively. opts provides confirmation of the 
% options utilised, together with the number of function 
% evaluations etc. (type help foptions for details). For the case of 
% AR(1) or SRW models, alpha less than unity specifies the 
% additional parameter, while the default value of unity implies a 
% RW or IRW model. Finally, parse are the standard errors of the 
% NVR and alpha (if optimised) hyper-parameters. However, 
% computation time can sometimes be greatly reduced if this 
% 4th output argument is omitted from the function call.

if nargin==0
  disp(' ')
  disp(' DLROPT  Hyper-parameter estimation for DLR')
  disp(' ')
  disp(' [nvr,alpha,opts,parse]=dlropt(y,z,TVP,meth,nvrc,alphac,nvr0,alpha0,opts,ALG,tab,P0)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, z=[]; end
if nargin<3, IRW=[]; end
if nargin<4, method=[]; end
if nargin<5, nvr=[]; end
if nargin<6, alpha=[]; end
if nargin<7, nvr0=[]; end
if nargin<8, alpha0=[]; end
if nargin<9, opts=[]; end
if nargin<10, ALG=[]; end
if nargin<11, output=[]; end
if nargin<12, P0=[]; end
if nargin<13, p=[]; end

if nargout>3
  [NVR,ALPHA,opts,separ]= dlropt0(y,z,IRW,method,nvr,alpha,nvr0,alpha0,opts,ALG,output,P0,p);
else
  [NVR,ALPHA,opts]= dlropt0(y,z,IRW,method,nvr,alpha,nvr0,alpha0,opts,ALG,output,P0,p);  
end

% end of m-file