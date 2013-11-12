% RIVDEMO1  Captain Toolbox demonstration
%
% Estimation of Transfer Function models from
% simulated input-output data
%
% See also RIV, RIVID, RIVDEMO2, RIVDEMO3

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% RIVDEMO1  Captain Toolbox demonstration
 
% This script analyses the output from a discrete time
% Transfer Function (TF) using the Refined Instrumental
% Variable (RIV) algorithm.
 
% We will start with a 2nd order TF with two samples time
% delay, one numerator parameter and a 2nd order noise model.
 
%                      b1(k)
%  y(k)  =  ---------------------------  u(k-2)
%            1 + a1.z^(-1) + a2.z^(-2) 
 
%                       1
%        +  ---------------------------  e(k)
%            1 + c1.z^(-1) + c2.z^(-2)
 
% where z^(-1) represents the backward shift operator,
% y(k) the output, u(k) the input and e(k) a zero mean
% white noise signal.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Transfer Function polynomials.   
 
a=[1 -1.75 0.8];  % denominator
b=[0 0 0.05];     % numerator
c=[1 -0.3 0.02];  % noise
 
% Simulation output.
 
u=[zeros(50, 1); ones(100, 1); zeros(50, 1)];  % input
en=[zeros(10, 1); filter(1, c, randn(190, 1)*0.05)];  % noise
y=filter(b, a, u) +  en;  % output
 
% Input-Output matrix
 
z=[y u]; plot(z)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Model structure.
 
na=2;  % denominator order
nb=1;  % numerator order
nd=2;  % pure time delay
nc=2;  % noise model
nn=[na nb nd nc];
 
% The Refined Instrumental Variable (RIV) algorithm is
% employed by using the following settings for 'flags'.
 
Ni=3;     % number of IV iterations (Least Squares if Ni=1)
Ft=1;     % filtering turned on
Nr=3;     % number of RIV iterations
flags=[Ni Ft Nr];  % see 'help riv' for details
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING MODEL : PLEASE WAIT
 
[TH, STATS, E]=riv(z, nn, flags);
 
% Plot results
 
subplot(211)
plot([u y y-E])
title('Input and Model fit')
subplot(212);
plot([en E]);
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
 
[b; bb; 0 0 sb(3)]  % numerator
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
[c; cc; 0 sb(4:5)]  % noise filter
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We will now try a multi-input example.
 
a=[1 -1.75 0.8];  % denominator
b1=[0 0 0.05];    % numerator 1
b2=[0 -0.1];      % numerator 2
c=[1 -0.6];       % noise
 
% Simulation output.
 
u1=[zeros(100, 1); ones(100, 1); zeros(50, 1)];  % input 1
u2=[zeros(50, 1); ones(100, 1); zeros(100, 1)];  % input 2
en=[zeros(10, 1); filter(1, c, randn(240, 1)*0.05)];  % noise
y=filter(b1, a, u1) + filter(b2, a, u2) + en;  % output
 
% Input-Output matrix
 
z=[y u1 u2]; clf; plot(z)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Model structure.
 
na=2;   % 2nd order denominator
nb1=1;  % numerator for input 1
nb2=1;  % numerator for input 2
nd1=2;  % pure time delay for input 1
nd2=1;  % pure time delay for input 2
nc=1;   % 1st order noise model
nn=[na nb1 nb2 nd1 nd2 nc];
 
% The Refined Instrumental Variable (RIV) algorithm is
% employed by using the following settings for 'flags'.
 
Ni=3;      % number of IV iterations (Least Squares if Ni=1)
Ft=1;      % filtering turned on
Nr=3;      % number of RIV iterations
Lr=1e-12;  % Linear solver switch
flags=[Ni Ft Nr Lr];  % see 'help riv' for details
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING MODEL : PLEASE WAIT
 
[TH, STATS, E]=riv(z, nn, flags);
 
% Plot results
 
subplot(211)
plot([u1 u2 y y-E])
title('Inputs and model fit')
subplot(212)
plot([E en])
title('''Actual'' and estimated observation noise')
 
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
 
[b1; bb(1, :); 0 0 sb(3)]  % numerator 1
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
[b2 0; bb(2, :); 0 sb(4) 0]  % numerator 2
 
[c; cc; 0 sb(5)]  % noise filter
 
echo off

% end of m-file