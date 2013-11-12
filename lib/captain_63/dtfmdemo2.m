% DTFMDEMO2  Captain Toolbox demonstration
%
% Dynamic Transfer Function (DTF) modelling
% of simulationed multi-input data
%
% See also DTFM, DTFMOPT, DTFMDEMO1, DARXDEMO

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for DTFMDEMO2: Paul McKenna

clear all
close all
format compact
echo on

clc
% DTFMDEMO2  Captain Toolbox demonstration
 
% This script analyses the output from a multi-input discrete 
% time Transfer Function (TF) using the functions for Dynamic
% TF modelling (DTF). Here, the model is not limited to the 
% ARX form required for DARX analysis.
 
% We will use two first order TFs each with the same gradually 
% changing time constant, one with a fixed numerator and a time 
% delay of three samples and the other numerator value switching 
% between two values and with a time delay of one sample.
 
%                     b1
%      y(k) =  ---------------- u1(k-3)
%              1 + a1(k).z^(-1)
 
%                      b2(k)
%              +  ---------------- u2(k-1)  +  e(k)
%                 1 + a1(k).z^(-1)
 
% where z^(-1) represents the backward shift operator, y(k) the 
% output, u(k) the input and e(k) a zero mean white noise 
% signal.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For this example, we utilise one constant numerator parameter 
% (b1=0.5) and a switching value for b2 (between 0.5 and -0.5) 
% and allow the a1 denominator parameter to change slowly over 
% time as a sine wave. Noise has been added to the output to 
% represent measurement errors, while Not-a-Number values are
% added in order to test the ability of the algorithms to
% handle interpolation.
 
load sdtfm2.dat
y=sdtfm2(:, 1);   % output (transfer function response)
u1=sdtfm2(:, 2);  % first input (white noise)
u2=sdtfm2(:, 3);  % second input (white noise)
a1=sdtfm2(:, 4);  % a1 denominator parameter
b1=sdtfm2(:, 5);  % b1 numerator parameter
b2=sdtfm2(:, 6);  % b2 numerator parameter
y0=sdtfm2(:, 7);  % noise free output 
y(1000:1100)=NaN*ones(101,1); % simulate missing data
y0(1000:1100)=NaN*ones(101,1); % simulate missing data
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The various signals are plotted below.
 
clf
subplot(411); plot([y y0]); ylabel('output');
title(['noise to signal variance ratio = ' ...
  int2str(100*var(y0(~isnan(y))-y(~isnan(y)))/var(y(~isnan(y)))) '%'])
subplot(412); plot(a1);ylabel('a_1')
subplot(413); plot(b1);ylabel('b_1')
subplot(414); plot(b2);ylabel('b_2')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We first call DTFMOPT to optimise the Noise Variance
% Ratio (NVR) hyper-parameters.
 
% The model structure is represented by the vector [A B T],
% where A is number of denominator parameters, B is the 
% vector of numbers of b parameters (one for each input) and
% T is the vector of total time delays (one for each input).
 
nn=[1 1 1 3 1];
 
% Since it is constant, we will fix the NVR for the first
% numerator parameter at zero. This will speed up the
% optimisation for the purposes of this demonstration.
 
nvrc=[-2 0 -2];  % optimise the NVRs for a1 and b2 only
 
% An Integrated Random Walk (IRW) model is used for the
% time variable a1 and b2 parameters, while the second
% parameter is constant, so a RW is more appropriate.
 
TVP=[1 0 1];
  
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING HYPER-PARAMETERS : PLEASE WAIT
 
nvr=dtfmopt(y, [u1 u2], nn, TVP, [], nvrc)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The DTFM function utilises the optimsed hyper-parameters.
 
% ESTIMATING MODEL PARAMETERS : PLEASE WAIT
 
[tfs, fit, fitse, par, parse, e]=dtfm(y, [u1 u2], nn, TVP, nvr);
 
% The final estimate of the constant second numerator parameter:
 
par(length(par), 2)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The estimated denominator parameter and standard errors
% are compared with the actual values.
 
clf
subplot(311);plot(a1, 'r'); hold on
plot(par(:, 1),'b')
plot(par(:, 1)+parse(:, 1), ':b')
plot(par(:, 1)-parse(:, 1), ':b')
axis([0 1500 -1.5 1.5])
subplot(312);plot(b1, 'r'); hold on
plot(par(:, 2),'b')
plot(par(:, 2)+parse(:, 2), ':b')
plot(par(:, 2)-parse(:, 2), ':b')
axis([0 1500 0.4 0.6])
subplot(313);plot(b2, 'r'); hold on
plot(par(:, 3),'b')
plot(par(:, 3)+parse(:, 3), ':b')
plot(par(:, 3)-parse(:, 3), ':b')
axis([0 1500 -1.5 1.5])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For comparision, the biased parameters based on DARX analysis
% are estimated.
 
% ESTIMATING MODEL PARAMETERS : PLEASE WAIT
 
[tfs, fit, fitse, parx, parsex]=darx(y, [u1 u2], nn, TVP, nvr);
 
% The new estimates are shown in the plot.
 
subplot(311);
plot(parx(:, 1), 'g')
title('Actual(red), MDTF-IV (blue) and DARX-LS (green)')
subplot(312);
plot(parx(:, 2), 'g')
subplot(313);
plot(parx(:, 3), 'g')
 
echo off

% end of m-file
