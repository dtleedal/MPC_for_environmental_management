% DLRDEMO  Captain Toolbox demonstration
%
% Dynamic Linear Regression (DLR) analysis
% of bivariate advertising data
%
% See also DLR, DLROPT

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% DLRDEMO  Captain Toolbox demonstration
 
% This script analyses advertising data using
% the functions for Dynamic Linear Regression (DLR).
 
load adv.dat
x=adv(:, 1);  % scaled expenditure on advertising
y=adv(:, 2);  % measure of response to advertising (0-1)
 
% We will use the following model: y(t) = c1 + c2(t)*x(t)
% The corresponding regressors are formed below.
 
z=[ones(size(x)) x];
 
% The response data contains missing values, represented
% in Matlab as Not-a-Number variables and forming gaps
% in the plot.
 
clf;
subplot(211); plot(y)
subplot(212); plot(x)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We first call DLROPT to optimise the Noise Variance
% Ratio (NVR) hyper-parameters.
 
% A Random Walk (RW) model is used for the
% time variable parameters.
 
TVP=0;
  
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING HYPER-PARAMETERS : PLEASE WAIT
 
nvr=dlropt(y, z, TVP)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The DLR function utilises the optimsed hyper-parameters.
 
[fit, fitse, par, parse]= dlr(y, z, TVP, nvr);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The parameter estimates and their standard errors
% are plotted.
 
clf
subplot(211)
plot(par(:, 1))
hold on
plot(par(:,1)+2*parse(:, 1),':')
plot(par(:,1)-2*parse(:, 1),':')
title('c1')
set(gca, 'xlim', [0 length(par)])
 
subplot(212)
plot(par(:, 2))
hold on
plot(par(:, 2)+2*parse(:, 2),':')
plot(par(:, 2)-2*parse(:, 2),':')
title('c2')
set(gca, 'xlim', [0 length(par)])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The data lie within the standard error bounds.
 
clf
plot(fit)
hold on
plot(fit+2*fitse, ':')
plot(fit-2*fitse, ':')
plot(y, 'o')
title('Data (o) and DLR model fit')
set(gca, 'xlim', [0 length(y)])
 
echo off

% end of m-file