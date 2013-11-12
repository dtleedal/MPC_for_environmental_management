% DARXDEMO  Captain Toolbox demonstration
%
% Dynamic Auto-Regression with eXogenous variables (DARX)
% modelling of simulated data
%
% See also DARX, DARXOPT, DTFMDEMO1, DTFMDEMO2

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for DARXDEMO: Paul McKenna

clear all
close all
format compact
echo on

clc
% DARXDEMO  Captain Toolbox demonstration
 
% This script analyses the output from a discrete time
% Transfer Function (TF) using the functions for Dynamic
% Auto-Regressive eXogenous variables (DARX) modelling.
 
% We will use a first order TF with a gradually increasing
% steady state gain, a fixed time constant and a time delay
% of two samples.
 
%              b(k)                       1
%  y(k) = -------------- u(k-2)  +  -------------- e(k)
%          1 + a.z^(-1)              1 + a.z^(-1)
 
% where z^(-1) represents the backward shift operator,
% y(k) the output, u(k) the input and e(k) a zero mean
% white noise signal.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For this example, we utilise a constant denominator
% parameter (a=-0.8) and allow the numerator parameter
% to ramp up slowly over time. Noise has been added to
% the output to represent measurement errors. The
% signals are loaded and plotted below.
 
load sdarx.dat
y=sdarx(:, 1);  % output (transfer function response)
u=sdarx(:, 2);  % input (white noise)
b=sdarx(:, 3);  % numerator parameter for reference
 
clf;
subplot(211)
plot([y u])
subplot(212)
plot(b)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We first call DARXOPT to optimise the Noise Variance
% Ratio (NVR) hyper-parameters.
 
% The model structure is represented by the triad [A B T],
% where A and B are the number of denominator and numerator
% parameters respectively, while T is the total time delay.
 
nn=[1 1 2];
 
% A Random Walk (RW) model is used for the
% time variable parameters.
 
TVP=0;

% Since it is constant, we will fix the NVR for the
% denominator parameter at zero. This will speed up the
% optimisation for the purposes of this demonstration.
 
nvrc=[0 -2];  % optimise the second NVR only
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING HYPER-PARAMETERS : PLEASE WAIT
 
nvr=darxopt(y, u, nn, TVP, [], nvrc)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The DARX function utilises the optimsed hyper-parameters.
 
[tfs, fit, fitse, par, parse]= darx(y, u, nn, TVP, nvr);
 
% The final estimate of the denominator parameter
% compares favourably with the actual value (a=-0.8).
 
par(length(par), 1)
 
% The estimated numerator parameter and standard errors are
% compared with the actual value (b increases from 4 to 5).
 
clf
plot([par(5:end, 2) b(5:end)])
hold on
plot(par(5:end, 2)+parse(5:end, 2), ':')
plot(par(5:end, 2)-parse(5:end, 2), ':')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% It is important to note that DARX provides Least Squares
% estimates of the parameters and suffers from the same
% limitations as conventional en-bloc algorithms, namely
% parameter bias when the noise does not conform to the ARX
% assumption; see e.g. Young, P.C. (1984) Recursive Estimation
% and Time Series Analysis, Springer-Verlag.
 
% In such cases, more sophisticated bivariate modelling tools
% are required, such as the Instrumental Variable aglorithms
% included in the following functions: DTFM, RIV and RIVID.
 
echo off

% end of m-file