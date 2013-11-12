function [amp, t, peaks, amps] = arspec(x, nar, opt, N, vec)
% ARSPEC  Auto-Regression spectrum
%
% [amp,t,pks,amps]=arspec(y,nar,out,N,v)
%
% y: Time series (*)
% nar: AR-spectrum order (0)
%        0: automatic identification based on AIC
%        Scalar: AR-spectrum order
%        Vector: specify AR polynomial
% out: Table and graphical output ([1 0]) 
%        out(1): graph and text on (1), no graphics (-1), both off (0)
%        out(2): log scale (0: off; other: maximum period allowed)
% N: Normalised frequency axis (in the range 0-0.5) or number of
%      frequencies ordinates to use in the estimation (1032)
% v: Vector of periods where vertical lines should be plotted
%
% amp: Estimated spectrum
% t: Frequency axis at which the spectrum is estimated
% pks: Periods at which peaks occur
% amps: Amplitude of the peaks
%
% See also AIC, MAR, PERIOD, DHROPT

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is the only compulsory input 
% to this function.
% 
% The second input argument nar is the desired AR model order 
% on which the spectrum estimation is based. If it is not supplied, 
% is an empty matrix or is zero, the AR order is automatically 
% selected via the AIC criterium (see aic). The next output 
% argument out, is a vector of dimension 2 that controls the 
% tabular and graphical output of the function, allows for a log 
% scale, or selects a range of frequencies on which the AR-
% spectrum should be estimated. In particular, a value of 1 for 
% out(1) sets the graphical and text output on, if it is -1 there is no 
% graphical output, while both graphics and text are off for a 0 
% value. The scale for the power axis is set always to logarithmical 
% but the frequency axis scale may be selected by setting out(2) to 
% 1. Any other value of out(2) is regarded as a desired period, 
% with the graphical output then limited between the range   and 
% out(2).
% 
% N may be a scalar or a vector. When scalar it refers to the 
% number of points in which the frequency axis should be divided 
% (default value is 1032). Alternatively, a vector would be 
% considered as the frequency axis itself, normalised between 0 
% and 0.5 (corresponding to periods from   to 2). In order to 
% make the visual inspection of the graphical output more 
% straightforward, v can be supplied as a vector of any size, 
% indicating the periods at which vertical lines should be plotted.
% 
% The function returns amp, the estimated AR spectrum; t, the 
% frequency axis at which the spectrum is estimated; pks, the 
% frequency at which peaks occur; and amps, the amplitude of the 
% peaks. All these outputs are shown by the optional tabular 
% output.

if nargin==0
  disp(' ')
  disp(' ARSPEC  Auto-Regression spectrum')
  disp(' ')
  disp(' [amp,t,pks,amps]=arspec(y,nar,out,N,v)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, nar=[]; end
if nargin<3, opt=[]; end
if nargin<4, N=[]; end
if nargin<5, vec=[]; end

[amp, t, peaks, amps]= arspec0(x, nar, opt, N, vec);

% end of m-file