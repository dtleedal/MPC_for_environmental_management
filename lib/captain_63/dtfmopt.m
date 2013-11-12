function [NVR,opts,separ] = dtfmopt(y,u,nn,IRW,method,nvr,nvr0,opts,ALG,output,P0)
% DTFMOPT  Hyper-parameter estimation for DTFM
%
% [nvr,opts,parse]=dtfmopt(y,u,nn,TVP,meth,nvrc,nvr0,opts,ALG,tab,P0)
%                          1 2 3   4   5    6    7    8    9  10  11
%
% y: Time series (*)
% u: Input (*)
% nn: Model structure [na,nb(1:nu),nd(1:nu)] ([1 1 0])
% TVP: Model type for each TVP (0-RW, 1-IRW) (0)
% meth: Estimation method ('ml')
%     'ml': Maximum Likelihood
%     'f#': Sum of squares of the #-step-ahead forecasting errors
% nvrc: Constraints for each NVR (-2)
%       -2: Free estimation
%       -1: Constrained estimation (all parameters with nvrc=-1 are equal)
%      >=0: NVR constrained to this value (it is not estimated)
% nvr0: Initial NVR hyper-parameters (0.0001)
% opts: Optimisation options, see help for OPTIMSET (or FOPTIONS for ALG=3)
% ALG: Optimisation algorithm: 0=fminsearch, 1=fminunc, 2=lsqnonlin (not ML) (0)
%        ALG=3 uses older FMINS based optimisation with interrupt button
% tab: Display: 0=none, 1=tabulate results, 2=update window (2)
% P0: Initial P matrix (1e5)
%
% nvr: Estimated NVR hyper-parameters
% opts: Returned optimisation options. Type 'help foptions' for details
% parse: Standard Error of hyper-parameters (omit to reduce computation time)
%
% Example: dtfmopt(y, [u1 u2], [1 1 1 2 3], 0, [], [0 -2])
%   transfer function type model y(k) = a1(k)*y(k-1) + b1(k)*u1(k-2) + b2(k)*u2(k-3)
%   with a RW model for all parameters and a1 assumed constant (NVR fixed at zero)
%
% See also DTFM, FCAST, STAND, FOPTIONS, FMINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for DTFMOPT: Paul McKenna

% The time series vector y (column) and input matrix u are 
% specified by the user. Each column of u represents an input 
% signal. The function automatically handles missing values in y. 
% 
% In fact, y may be appended with additional NaNs to forecast or 
% backcast beyond the original series, as long as appropriate 
% values for u are also specified. The remaining input arguments 
% are optional. The DTF model structure is defined by nn, which 
% takes the form [n m d] where, in transfer function terms, n and 
% m are the number of denominator and numerator parameters 
% respectively, while d is the number of samples time delay. A 
% first order model with unity time delay and one numerator 
% parameter [1, 1, 1] is utilised by default.
% 
% TVP is a vector specifying the model associated with each DTF 
% parameter, listed in order of each denominator parameter and 
% then the numerator parameters for each input.
% Choices include a RW model by default (0) 
% or a IRW model (1). meth is the estimation method, where the 
% default Maximum Likelihood 'ml' may be replaced by 'f#' to 
% compute the sum of squares of the #-step-ahead forecasting 
% errors. nvrc defines the constraints for each NVR, where -2 
% implies free estimation, -1 constrained estimation (all 
% parameters with nvrc = -1 are equal) and >=0 implies the 
% associated NVR is constrained to this value (it is not estimated). 
% Initial NVR hyper-parameters may be specified using nvr0 
% respectively.
% 
% Optimisation options may be set using opts (type help foptions 
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
% are less than the total number of parameters, then they are 
% automatically expanded to the correct dimensions by using the 
% final element of the specified input vector. For example, if the 
% DTF model has 3 parameters but TVP is defined as [1 0], then 
% TVP is automatically expanded to [1 0 0]. Similarly, a scalar P0 
% implies an identity matrix scaled by this value.
% 
% The function returns nvr, the vector of NVR hyper-parameters. 
% opts provides confirmation of the options utilised, together with 
% the number of function evaluations etc. (type help foptions for 
% details). Finally, parse are the standard errors of the NVR 
% hyper-parameters. However, computation time can sometimes 
% be greatly reduced if this 3rd output argument is omitted from 
% the function call.

if nargin==0
  disp(' ')
  disp(' DTFMOPT  Hyper-parameter estimation for DTFM')
  disp(' ')
  disp(' [nvr,opts,parse]=dtfmopt(y,u,nn,TVP,meth,nvrc,nvr0,opts,ALG,tab,P0)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, u=[]; end
if nargin<3, nn=[]; end
if nargin<4, IRW=[]; end
if nargin<5, method=[]; end
if nargin<6, nvr=[]; end
if nargin<7, nvr0=[]; end
if nargin<8, opts=[]; end
if nargin<9, ALG=[]; end
if nargin<10, output=[]; end
if nargin<11, P0=[]; end
separ=[];
if nargout>2
  [NVR,opts]= dtfmopt0(y,u,nn,IRW,method,nvr,nvr0,opts,ALG,output,P0);
else
  [NVR,opts,separ]= dtfmopt0(y,u,nn,IRW,method,nvr,nvr0,opts,ALG,output,P0);
end

% end of m-file