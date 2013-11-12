% SDPDEMO  Captain Toolbox demonstration
%
% State Dependent Parameter (SDP) analysis of simulated data
%
% See also SDP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for SDPDEMO: Paul McKenna

clear all
close all
format compact
echo on

clc
% SDPDEMO  Captain Toolbox demonstration
 
% This script demonstrates the use of the State Dependent
% Parameter (SDP) routine for obtaining a non-parametric
% state dependent model from a simulated dataset.
 
% For this example, the output is formed from the sum of
% three inputs (u1, u2 and u3), each multiplied by its own 
% state dependent parameter (a, b and c).
 
% y(k) = a(k) * u1(k)  +  b(k) * u2(k)  +  c(k) * u3(k)  
 
% The SDPs themselves are each dependent on the given states 
% (x1, x2 and x3). In the following example, they are simply
% white noise signals with different seeds. If no dependent
% states are supplied by the user, the algorithm assumes that
% each parameter is dependent upon its associated input regressor,
% (so that in this example a would be dependent on u1, b on u2).
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% length of the simulated data set
 
n=1500;
 
% dependent states
 
randn('seed',0); x1=randn(n,1);  % for first regressor
randn('seed',1); x2=randn(n,1);  % for second regressor
randn('seed',2); x3=randn(n,1);  % for third regressor
x=[x1 x2 x3];
 
% state dependent parameters
 
a=(0.5*x1)+3;       % linear state dependency for the first SDP
b=(x2.^2)+1;        % quadratic state dependency for the second SDP
c=(sin(x3*pi/2))+2; % sinusoidal state dependency for the third SDP
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The inputs signals are more white noise sequences.
 
randn('seed',3); u1=randn(n,1);  % first input
randn('seed',4); u2=randn(n,1);  % second input
randn('seed',5); u3=randn(n,1);  % third input
z=[u1 u2 u3];
 
% The output is then generated from the inputs and their
% associated state dependent parameters.
 
y=(a.*u1)+(b.*u2)+(c.*u3);  % output
 
% Missing values are added in order to test the ability of the
% algorithms to automatically handle interpolation
 
u1(200:220)=NaN;
y(400:420)=NaN;
 
% Finally, white noise is added to the output signal to
% represent measurement noise.
 
randn('seed',6); ng=randn(n,1); ng=3*ng;
y=y+ng;
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The noise to signal ratio is high, as seen by plotting the 
% uncontaminated data and the noise signal together.
 
subplot(211)
plot(y-ng,'-k'); ax=axis; ylabel('Output (pre-noise)')
 
title([' Noise to signal ratio = ' ...
  int2str(100*covnan(ng)/covnan(y-ng)) '%  by variance']);
 
subplot(212)
plot(ng,'-k'); axis(ax); ylabel('Noise signal')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The output, inputs and states are all plotted.
 
clf
subplot(4,2,1); plot(y,'-r');ylabel('Output')
subplot(4,2,3); plot(u1,'-b');ylabel('Input 1')
subplot(4,2,5); plot(u2,'-b');ylabel('Input 2')
subplot(4,2,7); plot(u3,'-b');ylabel('Input 3')
subplot(3,2,2); h=plot(x1,a,'.m');
xlabel('Dependent state 1'); ylabel('SDP 1'); set(h,'markersize',1)
subplot(3,2,4); h=plot(x2,b,'.m');
xlabel('Dependent state 2'); ylabel('SDP 2'); set(h,'markersize',1)
subplot(3,2,6); h=plot(x3,c,'.m');
xlabel('Dependent state 3');ylabel('SDP 3');set(h,'markersize',1)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The NVRs for the SDPs can be provided explicitly or 
% optimised by an internal call to the DLROPT function from
% SDP. If the any of the given NVRs are negative, this is 
% taken to indicate that the user wishes to optimise the NVR
% for the corresponding SDP. An NVR of -1 will force the 
% algorithm to optimise at the first iteration for that SDP,
% then retain this value for subsequent iterations. In 
% general an NVR input of -n will cause optimisation at the 
% first n iterations.
 
% In the following example, we optimise the first NVR at the
% first two iterations and the other NVRs at only the first
% iteration.
 
nvr=[-2 -1];  
 
% A Random Walk (RW) model is used for the first SDP and
% an Integrated Random Walk (IRW) model for all other SDPs.
 
TVP=[0 1];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The estimation options are set such that there can be a
% maximum of 10 iterations, the convergence criteria (for
% r-squared) is 0.0001, the plotting option is set to 1 
% to show results during estimation and the other estimation
% options are set to their default (shown by -1 inputs).
 
opts=[10 0.0001 -1 -1 -1 1];
 
% The main function call is then made, to carry out
% the estimation, plotting the estimates as they are made.
% The initial esimates with NVR=0 are plotted blue. The
% latest SDP estimates are plotted red, with the points from
% the previous iterations being grey.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
clf
[fit,fitse,par,parse,zs,pars,parses,rsq,nvrid,y0]=sdp(y,z,x,TVP,nvr,opts);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Identified NVR for the first SDP
 
nvrid(1)
 
% Identified NVR for the second SDP
 
nvrid(2)
 
% Identified NVR for the third SDP
 
nvrid(3)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The r-squared value between the model output and the 
% output data is returned.
 
rsq
 
% Finally, the results can be plotted against the true SDP
% relationships, with standard error bounds shown as dashed
% lines.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
figure(1); clf
l1=plot(x1,a,'.b',zs(:,1),pars(:,1),'-g'); hold on
set(l1(2),'linewidth',2); set(l1(1),'markersize',1)
ylabel('a_k'); xlabel('x1_k');axis('square'); legend('par','sdp',-1)
plot(zs(:,1),pars(:,1)+parses(:,1),'--g',zs(:,1),pars(:,1)-parses(:,1),'--g')
title('First state dependent parameter')
 
figure(2)
l1=plot(x2,b,'.b',zs(:,2),pars(:,2),'-g'); hold on
set(l1(2),'linewidth',2); set(l1(1),'markersize',1)
ylabel('b_k'); xlabel('x2_k'); axis('square'); legend('par','sdp',-1)
plot(zs(:,2),pars(:,2)+parses(:,2),'--g',zs(:,2),pars(:,2)-parses(:,2),'--g')
title('Second state dependent parameter')
 
figure(3)
l1=plot(x3,c,'.b',zs(:,3),pars(:,3),'-g'); hold on
set(l1(2),'linewidth',2);set(l1(1),'markersize',1)
ylabel('c_k'); xlabel('x3_k'); axis('square'); legend('par','sdp',-1)
plot(zs(:,3),pars(:,3)+parses(:,3),'--g',zs(:,3),pars(:,3)-parses(:,3),'--g')
title('Third state dependent parameter')
 
echo off

% end of m-file