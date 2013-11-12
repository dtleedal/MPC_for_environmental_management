% DARDEMO  Captain Toolbox demonstration
%
% Dynamic Auto-Regression (DAR) analysis of a
% signal with sawtooth changing frequency
%
% See also DAR, DAROPT, DARSP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% DARDEMO  Captain Toolbox demonstration
 
% This script analyses a signal with sawtooth
% changing frequency using the functions for Dynamic
% Auto-Regression (DAR). The signal is loaded below.
 
load sdar.dat

% Data should be standardised before DAR analysis.
 
y=stand(sdar);
clf; plot(y);

% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% We first call DAROPT to optimise the Noise Variance
% Ratio (NVR) hyper-parameters.
 
% We will fit a 2nd order autoregression model.
 
nar=[1:2];

% The first parameter will be modelled with an Integrated
% Random Walk (IRW) and the second with a RW.
 
TVP=[1 0];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING HYPER-PARAMETERS : PLEASE WAIT
 
nvr=daropt(y, nar, TVP)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% DAR returns the model components and parameters.
 
[fit, fitse, par, parse]= dar(y, nar, TVP, nvr);
 
% The parameters and their standard errors are shown.
 
clf; plot(par)
hold on
plot(par+2*parse, ':')
plot(par-2*parse, ':')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Finally, the time varying AR spectrum is plotted.
 
clf
darsp(par, 6, 1);
 
% The 2nd and 3rd input arguments specify that a 3D spectrum
% plot with a resolution of 2^6 is required. The visual
% appearance of this plot depends on the computer platform
% and it is sometimes necessary to experiment to find the
% best resolution. Contoured surfaces and 2D stacked plots
% may be obtained by changing the 3rd input argument.
 
echo off

% end of m-file