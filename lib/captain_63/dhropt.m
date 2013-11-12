function [bl,alpha,opts,amp,separ] = dhropt(y,P,IRWharm,nar,nvr,alpha,nvr0,alpha0,opts,ALG,output,t,Interv)
% DHROPT  Hyper-parameter estimation for DHR
%
% [nvr,alpha,opts,amp,parse]=dhropt(y,P,TVP,meth,nvrc,alphac,nvr0,alpha0,opts,ALG,tab,tf,Int)
%                                   1 2  3   4    5     6     7     8     9   10  11  12 13
%
% y: Time series (*)
% P: Periodic components; set P(1)=0 to include a trend (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW, 2-Trigonommetric) (0)
%        (for LLT use RW and IRW trends simultaneously)
% meth: Estimation method (24)
%         meth>0: Frequency domain optimisation based on the AR(meth) spectrum
%         meth<0: ML in time domain with abs(meth) used as the order of the AR 
%                 spectrum to estimate initial conditions for hyper-parameters
%         meth='f#': Sum of squares of the #-steps-ahead forecasting errors with
%                    initial conditions from the frequency domain using AR(24)
%         meth='f# n': As above but specifying AR(n) to estimate initial conditions
% nvrc: Constraints for each NVR (-2)
%     -2: Free estimation
%     -1: Constrained estimation (all parameters with nvrc=-1 are equal)
%    >=0: NVR constrained to this value (it is not estimated)
% alphac: Constraints for each alpha (-2, -1, or >=0 as for nvrc) (1)
% nvr0: Initial NVR hyper-parameters (0)
% alpha0: Initial alpha hyper-parameters (1)
% opts: Optimisation options, see help for OPTIMSET (or FOPTIONS for ALG=3)
% ALG: Optimisation algorithm: 0=fminsearch, 1=fminunc, 2=lsqnonlin (not ML) (2)
%        ALG=3 uses older FMINS based optimisation with interrupt button
% tab: Display: 0=none, 1=tabulate results, 2=update window and frequency plot (2)
% tf: User supplied frequency information (1032)
%       scalar: Number of subdivisions of frequency axes
%       vector: Frequency axis (between 0-0.5)
%       matrix (nx2): Frequency axis (1st column) and AR spectrum amplitude
% Int: Vector of variance intervention points (0)
%
% nvr: Estimated NVR hyper-parameters
% alpha: Estimated alpha hyper-parameters
% opts: Returned optimisation options. Type 'help foptions' for details
% amp: Spectra [t, amp, ampm], e.g. use semilogy(t, amp, t, ampm)
%        t: Frequency axis at which the spectra are evaluated
%        amp: Empirical spectrum
%        ampm: Fitted model spectrum
% parse: Standard Error of hyper-parameters (omit to reduce computation time)
%
% Example: dhr(y, [0 12./(1:6)], [1 0], 24, -2, [-2 1])
%   optimise for a SRW trend, together with 6 periodic components
%   (12 and harmonics) each modelled with a RW (alpha fixed at 1)
%
% See also DHR, FCAST, STAND, FOPTIONS, FMINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) and associated m periodic 
% components P are specified by the user. Set the first element of 
% P to include a trend. For example, [0 12] implies a trend and a 
% seasonal component for monthly data. The function 
% automatically handles missing values in y. In fact, y may be 
% appended with additional NaNs to forecast or backcast beyond 
% the original series. The remaining input arguments are optional.
% 
% TVP is a vector specifying the model associated with each 
% regression parameter, listed in the same order as the elements of 
% P. Choices include a RW/AR(1) model by default (0) or a 
% IRW/SRW model (1). The 4th input argument meth selects the 
% estimation method, where the default frequency domain 
% optimisation based on the AR(24) spectrum may be replaced by 
% an AR(meth) spectrum by specifying a positive scalar. A 
% negative scalar implies Maximum Likelihood in time domain, 
% with -meth used as the order of the AR spectrum to estimate the 
% initial conditions. Finally, specifying 'f#' computes the sum of 
% squares of the #-step-ahead forecasting errors, with the initial 
% conditions obtained from a frequency domain optimisation using 
% the AR(24) spectrum, or 'f# n' to specify the AR(n) spectrum.
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
% 
% When meth = 'f#' there is the additional option of using leastsq 
% (see help leastsq). Note that if ALG = 2 and meth = 'ml', then 
% an error occurs, since leastsq cannot be used in the Maximum 
% Likelihood case. The Optimisation Toolbox for MATLAB(r) is 
% required to use fminu or leastsq.
% 
% To set the display options, if tab = 1 or 2, then the final results 
% are displayed in tabular form. Additionally, if tab = 2, a window 
% appears during optimisation showing the latest value of the 
% Likelihood Function or the Sum-of-Squares for the #-step-ahead 
% forecasting errors. When ALG = 0 (fmins) and tab = 2, a 'stop 
% button' will appear below the update window: click to terminate 
% the optimisation and return the current estimates.
% 
% For frequency domain optimisation, tf may be a scalar or a 
% vector. When scalar it refers to the number of points in which 
% the frequency axis should be divided (default value is 1032). 
% Alternatively, a vector would be considered as the frequency 
% axis itself, normalised between 0 and 0.5 (corresponding to 
% periods from inf to 2). Alternatively, a matrix specifies the 
% frequency axis (1st column) and AR spectrum amplitude 
% directly.
% 
% Finally, Int allows for sharp (discontinuous) local changes in the 
% parameters at the user supplied intervention points. These need 
% to be defined either manually or by some detection method for 
% sharp local changes. Here, Int should take the same dimensions 
% as y, with positive values indicating variance intervention 
% required. Note that this option does not influence the NVR 
% estimates when meth>0 (the default option).
% If the lengths of TVP, nvrc, alphac, nvr0, alpha0 or P0 or x0 
% are less than m, then they are automatically expanded to the 
% correct dimensions by using the final element of the specified 
% input vector. For example, if z has 3 columns but TVP is 
% defined as [1 0], then TVP is automatically expanded to [1 0 0]. 
% Similarly, a scalar P0 implies an identity matrix scaled by this 
% value.
% 
% The function returns vectors of NVR and   hyperparameters, 
% nvr and alpha respectively. opts provides confirmation of the 
% options utilised, together with the number of function 
% evaluations etc. (type help foptions for details). For the case of 
% AR(1) or SRW models, alpha less than unity specifies the 
% additional parameter, while the default value of unity implies a 
% RW or IRW model. amp returns the spectra used in the 
% optmisation in the form [t, am, ampm], where t  is the 
% frequency axis at which the spectra are evaluated, am is the 
% empirical spectrum and ampm is the fitted model spectrum. The 
% associated spectrum can be graphed with e.g. semilogy(t, amp, 
% t, ampm).
% 
% Finally, parse are the standard errors of the NVR and alpha (if 
% optimised) hyper-parameters. However, computation time can 
% sometimes be greatly reduced if this 4th output argument is 
% omitted from the function call.

if nargin==0
  disp(' ')
  disp(' DHROPT  Hyper-parameter estimation for DHR')
  disp(' ')
  disp(' [nvr,alpha,opts,amp,parse]=dhropt(y,P,TVP,meth,nvrc,alphac,nvr0,alpha0,opts,ALG,tab,tf,Int)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, P=[]; end
if nargin<3, IRWharm=[]; end
if nargin<4, nar=[]; end
if nargin<5, nvr=[]; end
if nargin<6, alpha=[]; end
if nargin<7, nvr0=[]; end
if nargin<8, alpha0=[]; end
if nargin<9, opts=[]; end
if nargin<10, ALG=[]; end
if nargin<11, output=[]; end
if nargin<12, t=[]; end
if nargin<13, Interv=[]; end

if nargout>4
  [bl,alpha,opts,amp,separ]=dhropt0(y,P,IRWharm,nar,nvr,alpha,nvr0,alpha0,opts,ALG,output,t,Interv);
else
  [bl,alpha,opts,amp]=dhropt0(y,P,IRWharm,nar,nvr,alpha,nvr0,alpha0,opts,ALG,output,t,Interv);
end

% end of m-file