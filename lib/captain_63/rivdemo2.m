% RIVDEMO2  Captain Toolbox demonstration
%
% Identification of Transfer Function models from
% simulated input-output data
%
% See also RIV, RIVID, RIVDEMO1, RIVDEMO3

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% RIVDEMO2  Captain Toolbox demonstration
 
% This script analyses the output from a discrete time
% Transfer Function (TF) using the Refined Instrumental
% Variable (RIV) algorithm.
 
% We will start with a 2nd order TF with three samples
% time delay and one numerator parameter.
 
%                      b1(k)
%  y(k)  =  ---------------------------  u(k-3)  +  e(k)
%            1 + a1.z^(-1) + a2.z^(-2) 
 
% where z^(-1) represents the backward shift operator,
% y(k) the output, u(k) the input and e(k) a zero mean
% white noise signal.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Transfer Function polynomials.   
 
a=[1 -1.75 0.8];  % denominator
b=[0 0 0 0.05];  % numerator
 
% Simulation output.
 
u=[zeros(50, 1); ones(100, 1); zeros(50, 1)];  % input
e=[zeros(10, 1); randn(190, 1)*0.005];  % white noise
y=filter(b, a, u) + e;  % output
 
% Input-Output matrix
z=[y u];
plot(z)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We will use RIVID to determine a satisfactory model
% structure, searching for up to 3 parameters each for
% the numerator and denominator polynomials and for a
% possible 0 to 3 time delays. There is no noise model.
 
nn=[1 1 0 0; 3 3 3 0]
 
% The Refined Instrumental Variable (RIV) algorithm is
% employed by using the following settings for 'flags'.
 
Ni=3;     % number of IV iterations (Least Squares if Ni=1)
Ft=2;     % filtering turned on
Nr=0;     % number of RIV iterations
flags=[Ni Ft Nr];  % see 'help riv' for details
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING MODEL : PLEASE WAIT
 
[TH, STATS, E]=rivid(z, nn, [], flags);
 
% Scroll up to see the top 20 models in terms of their YIC
% value, where the lowest YIC at the top of the list indicates
% the best compromise between a good fit to the data and a
% high confidence in the parameter estimates (low variance).
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Plot results
 
figure(1); clf
subplot(211)
plot([u y y-E])
title('Input and Model fit')
subplot(212);
plot([e E]);
title('''Actual'' and estimated coloured observation noise')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Model fit, where unity implies a perfect fit.
 
RT2=STATS(3)  % Coefficient of Determination
 
% Extract parameter estimates.
 
[aa, bb, cc, P0]=getpar(TH);
 
% The covariance matrix P0 provides the standard errors.
 
sb=sqrt(diag(P0))';
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The parameter estimates are displayed below. In each case,
% the actual values are shown in the first row, followed by
% the estimates and finally the standard deviations.
 
[a; aa; 0 sb(1:2)]  % denominator
 
[b; bb; 0 0 0 sb(3)]  % numerator
 
echo off

% end of m-file