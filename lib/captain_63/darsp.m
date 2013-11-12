function [H,ph] = darsp(par, res, PL, Nc, dt)
% DARSP  Dynamic Auto-Regression spectra plot
%
% [H,ph]=darsp(par,res,PL,Nc,dt)
%
% par: Non-stationary parameters from DAR (*)
% res: Resolution of H in log base 2 (6 - equivalent to 2^6)
% PL: Plot type (0)
%       0: No plot (H and ph are returned)
%       1: 3d coloured surface with contours
%       2: 2d contour plot
%       3: Stacked plot
% Nc: Number of contours (PL=2) or stacking distance (PL=3) (10)
% dt: Sampling rate (e.g. seconds) (1)
%
% H: amplitude of AR spectrum (e.g. use surf(H) for presentation)
% ph: phase (in radians, multiply by 180/pi to obtain degrees)
%
% Example: [fit, comp, par]=dar(y); darsp(par, 1)
%
% See also DAR, DAROPT, AIC, ARSPEC

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The parameter matrix par from DAR is the only compulsory 
% input to this function.
% 
% Use PL to specify the type of time spectra graph. The default 
% zero calculates the spectra and returns the output arguments but 
% does not display the graph. A value of 1 graphs a 3d coloured 
% surface with contours, 2 a 2d contour plot and 3 a stacked plot. 
% 
% The associated input arguments res and Nc control the 
% dimension of H, namely 2^res by length(par), and the number 
% of contours (PL = 2) or stacking distance (PL = 3) respectively. 
% 
% Experiment with res and Nc to obtain the best plot for different 
% computer systems and data sets. The value of Nc is ignored 
% when PL = 0 or PL = 1. Finally, dt is the sampling rate (e.g. 
% seconds) used when plotting the graph.
% 
% The outputs arguments H and ph are the DAR spectrum (for 
% example, use surf(H) to plot the time spectra as a surface) and 
% phase (in radians) respectively.

if nargin==0
  disp(' ')
  disp(' DARSP  Dynamic Auto-Regression spectra plot')
  disp(' ')
  disp(' [H,ph]=darsp(par,res,PL,Nc,dt)')
  disp(' ')
  return
end

if nargin<1, par=[]; end
if nargin<2, res=[]; end
if nargin<3, PL=[]; end
if nargin<4, Nc=[]; end
if nargin<5, dt=[]; end

if nargout==2
   [H,ph]=darsp0(par, res, PL, Nc, dt);
elseif nargout==1
   [H]=darsp0(par, res, PL, Nc, dt);
elseif nargout==0
   darsp0(par, res, PL, Nc, dt);
end

% end of m-file