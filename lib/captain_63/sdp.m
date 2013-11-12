function [fit,fitse,par,parse,zs,pars,parses,rsq,nvr,y0] = sdp(y,z,x,TVP,nvr,opts,P0,x0,nvr0,tab)
% SDP  Non-parametric state dependent modelling (backfitting)
%
% [fit,fitse,par,parse,xs,pars,parses,rsq,nvre,y0]=sdp(y,z,x,TVP,nvr,opts,P0,x0,nvr0,tab)
%                                                      1 2 3  4   5   6   7   8   9  10
%
% y: Time series (*)
% z: Regressors (*)
% x: States on which the parameters depend, column for each regressor (z)
% TVP: Model type for each TVP (0-RW, 1-IRW) (0)
% nvr: NVR hyper-parameters (0)
%        if an nvr is a negative integer (-n) it is
%        optimised over the first n backfit steps
% opts: estimation options [iter con meth sm ALG plotopt], use -1 for defaults
%         iter: number of backfitting iterations (10)
%         con: backfitting convergence threshold (0.001)
%         meth: Optimisation estimation method (0)
%                 0: Maximum Likelihood
%                 Integer:Sum of squares of the #-step-ahead forecasting errors
%         sm: Smoothing on (1-default) or off (0-saves memory)
%         ALG: Smoothing algorithm P (0) or Q (1-default)
%         plotopt: Plot results during estimation on (1) or off (0-default)
% P0: Initial P matrix diagonal (1e5)
% x0: Initial state vector (0)
% nvr0: Initial NVRs specified by user (0.0001)
% tab: Display: 0=none, 1=tabulate results, 2=update window (2)
%
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% xs: Sorted states
% pars: Sorted parameter estimates
% parses: Sorted standard errors of parameters
% rsq: R squared value
% nvre: Estimated NVRs
% y0: Interpolated data
%
% Example: sdp(y,[u1 u2],[x1 x2],[0 1],[-2 -1])
%   regression type model y = c1(x1)*u1 + c2(x1)*u2, with an RW model for
%   c1 where the dependent state is x1 and an IRW model for c2 where the 
%   dependent state is x2; the NVR for the first SDP (c1) is optimised at
%   the first iteration and at the first two iterations for the second SDP
%
% See also FCAST, STAND, DLR, DHR

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for SDP: Paul McKenna

% The time series vector y (column), together with the associated 
% m regressors z and states x are specified by the user. Here, z and 
% x have the same number of rows as y, with a column for each 
% regressor/state. The function automatically handles missing 
% values in y. In fact, y may be appended with additional NaNs to 
% forecast or backcast beyond the original series. The remaining 
% input arguments are optional.
% 
% TVP is a vector specifying the model associated with each 
% regression parameter, listed in the same order as the columns of 
% z. Choices are a RW model by default (0) or a IRW model (1). 
% 
% nvr is a vector of NVR hyperparameters for each regressor 
% where, for example, zero (default) implies time invariant 
% parameters. Negative values imply that the NVR 
% hyperparameters are automatically optimised over the first -nvr 
% backfitting steps. opts is a vector of options: [iter con meth sm 
% ALG plotopt].
% 
% Here, iter is the number of backfitting iterations (default 10); 
% con is the backfitting convergence threshold (0.001); meth is 
% the optimisation estimation method, which may be either ML (0 
% - default) or, for a positive integer, the sum of squares of the 
% meth-step-ahead forecasting errors; sm specifies whether FIS 
% smoothing is on (1 - default) or off (0 - here, the model fit and 
% estimated parameters are their filtered values, which speeds up 
% the algorithm and reduces memory usage in cases when 
% smoothing is not required); ALG selects either the P (0) or Q (1 
% - default) smoothing algorithm (should convergence problems 
% be encountered, changing the algorithm in this manner may 
% help); and, finally, plotopt selects graphical display of the 
% results during estimation (1) or turns this option off (0 - 
% default). If only some of the values in opts require changing, the 
% rest may be set to -1 for their default values.
% 
% The initial state vector and diagonal of the P-matrix may be 
% specified using x0 and P0, with default values of 0 and 1e5 
% respectively. Finally, nvr0 is a vector of initial NVR 
% hyperparameters, utilised by the first backfitting step (0.0001).
% 
% If the lengths of TVP, nvr, P0 or x0 are less than m, then they 
% are automatically expanded to the correct dimensions by using 
% the final element of the specified input vector. For example, if z 
% has 3 columns but TVP is defined as [1 0], then TVP is 
% automatically expanded to [1 0 0]. Similarly, a scalar P0 implies 
% an identity matrix scaled by this value.
% 
% The function returns the model fit (with the same dimensions as 
% y) and parameters par (one column for each regressor), together 
% with the associated standard errors in each case, fitse and parse. 
% It also returns the vectors of sorted states, parameters and 
% standard errors, xs, pars and parses respectively. Finally, rsq is 
% a measure of model fit R2, nvre are the final NVR estimates 
% and y0 the interpolated data, where the latter consist of the 
% original series with any missing data replaced by the model.

if nargin==0
  disp(' ')
  disp(' SDP  Non-parametric state dependent modelling (backfitting)')
  disp(' ')
  disp(' [fit,fitse,par,parse,xs,pars,parses,rsq,nvre,y0]=sdp(y,z,x,TVP,nvr,opts,P0,x0,nvr0)')
  disp(' ')
  return
end

if nargin<10, tab=[]; end
if nargin<9, nvr0=[]; end
if nargin<8, x0=[]; end
if nargin<7, P0=[]; end
if nargin<6, opts=[]; end 
if nargin<5, nvr=[]; end
if nargin<4, TVP=[]; end
if nargin<3, x=[]; end
if nargin<2, z=[]; end
if nargin<1, y=[]; end

[fit,fitse,par,parse,zs,pars,parses,rsq,nvr,y0]=sdp0(y,z,x,TVP,nvr,opts,P0,x0,nvr0,tab);

% end of m-file