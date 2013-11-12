% DTFMDEMO1  Captain Toolbox demonstration
%
% Dynamic Transfer Function (DTF) modelling
% of simulated input-output data
%
% See also DTFM, DTFMOPT, DTFMDEMO2, DARXDEMO

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for DTFMDEMO1: Paul McKenna

clear all
close all
format compact
echo on

clc
% DTFMDEMO1  Captain Toolbox demonstration
 
% This script analyses the output from a discrete time
% Transfer Function (TF) using the functions for Dynamic
% TF modelling (DTF). Here, the model is not limited to
% the ARX form required for DARX analysis.
 
% We will use a first order TF with a gradually changing
% time constant, a fixed numerator and a time delay of
% two samples.
 
%                     b
%      y(k) = ---------------- u(k-2)  +  e(k)
%              1 + a(k).z^(-1)
 
% where z^(-1) represents the backward shift operator,
% y(k) the output, u(k) the input and e(k) a zero mean
% white noise signal.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For this example, we utilise a constant numerator
% parameter (b=0.5) and allow the denominator parameter
% to change slowly over time as a sine wave. Noise has
% been added to the output to represent measurement errors.
% The signals are loaded and plotted below.
 
load sdtfm1.dat
y=sdtfm1(:, 1);  % output (transfer function response)
u=sdtfm1(:, 2);  % input (white noise)
a=sdtfm1(:, 3);  % denominator parameter for reference
 
clf
subplot(211); plot(y)
subplot(212); plot(a)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We first call DTFMOPT to optimise the Noise Variance
% Ratio (NVR) hyper-parameters.
 
% The model structure is represented by the triad [A B T],
% where A and B are the number of denominator and numerator
% parameters respectively, while T is the total time delay.
 
nn=[1 1 2];
 
% Since it is constant, we will fix the NVR for the
% numerator parameter at zero. This will speed up the
% optimisation for the purposes of this demonstration.
 
nvrc=[-2 0];  % optimise the first NVR only
 
% An Integrated Random Walk (IRW) model is used for the
% time variable denominator parameter, while the second
% parameter is constant, so a RW is more appropriate.
 
TVP=[1 0];
  
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING HYPER-PARAMETERS : PLEASE WAIT
 
nvr=dtfmopt(y, u, nn, TVP, [], nvrc)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The DTFM function utilises the optimsed hyper-parameters.
%
% ESTIMATING MODEL PARAMETERS : PLEASE WAIT
 
[tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, nn, TVP, nvr);
 
% The final estimate of the numerator parameter:
 
par(length(par), 2)
 
% The estimated denominator parameter and standard errors
% are compared with the actual values.
 
clf
plot(a, 'm')
hold on
plot(par(:, 1))
plot(par(:, 1)+parse(:, 1), ':')
plot(par(:, 1)-parse(:, 1), ':')
axis([0 1500 -1.5 1.5])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For comparision, the biased parameters based on DARX analysis
% are estimated.
 
% ESTIMATING MODEL PARAMETERS : PLEASE WAIT
 
[tfs, fit, fitse, parx, parsex]=darx(y, u, nn, TVP, nvr);
 
% The new estimates are shown in the plot.
 
plot(parx(:, 1), 'r')
plot(parx(:, 1)+parsex(:, 1), ':r')
plot(parx(:, 1)-parsex(:, 1), ':r')
title('Actual, DTF(IV) and DARX(LS) (red)')

echo off

% end of m-file
