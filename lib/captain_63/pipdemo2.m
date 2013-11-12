% PIPDEMO2  Captain Toolbox demonstration
%
% Univariate Proportional-Integral-Plus (PIP)
% control design for global CO2 emissions
%
% See also PIPDEMO1, PIPDEMO3

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% PIPDEMO2  Captain Toolbox demonstration
 
% This script develops a univariate Proportional-Integral-Plus (PIP)
% control system to control global CO2 emissions, in order to
% stabilize atmospheric CO2 concentration levels in the atmosphere,
% based on a continuous-time Transfer Function (TF) model of the
% emissions-atmospheric CO2 data.
 
% The Matlab Control System Toolbox is required for this demo.
% If this is not installed on your system, the script will crash
% (or press Ctrl+C now to quit).
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The model used is based on an analysis of annual globally
% averaged data for the emissions (arising from fossil fuel 
% burning and land-use change) and atmospheric CO2 concentration.
% These are first loaded into the workspace as follows.
 
load globalco2.dat
u=globalco2(:, 2);
y=globalco2(:, 1);
 
plot(u, 'k'); hold on; plot(y, 'b')
title('Global CO2 (blue) and emissions (black)')

% A continuous-time TF model will be identified using RIVC. See
% RIVCDEMO for a worked example of such analysis and 'help RIVC'
% for a description of the input arguments.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ** Please wait while model is estimated **
[th, stats, e]=rivc([y-y(1) u-u(1) ones(size(y))], [1 1 1 5 5]);
plot(y-e, 'linewidth', 2);
 
% The estimated model (shown as the thick trace on the figure) is:
%
%             beta(1)                    1
%   y(t) = ------------- u(t-5)  +  ------------- beta(2)  +  e(k)
%            s + alpha                s + alpha
%
% where s represents the differential operator, y(t) the output 
% (atmospheric pCO2 level), u(t) the input CO2 emissions and e(k)
% is the residual noise. The parameter beta(2) is a constant 
% arising from the second 'constant' input to allow for non-zero
% initial condtions. The estimates are:
 
[alpha, beta]=getpar(th)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% In order to design a digital PIP control system, the continuous-
% time model has to be converted to discrete-time using the C2DM
% tool in the Matlab Control Systems Toolbox, ignoring the 5 period
% time delay (C2DM cannot handle this). Note that, because of this,
% the RIVC tool returns the numerator polynomial without any indication
% of the time delay, which the user must remember. This is because
% the approach of using leading zero elements in the numerator
% polynomial over the time delay length, to account for a time delay
% in discrete-time models, does not apply for continuous-time models
% (where the time delay is in time units not sampling intervals).
 
[b, a]=c2dm(beta(1), alpha, 1, 'zoh')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Now the PIPOPT function is used to compute the PIP controller
% coefficients for the implementation of the controller in TF form.
% See PIPDEMO1 for more details. Note that the 5 sample delay is
% now incorporated because we are considering the discrete-time
% equivalent of the continuous-time model.
 
v=pipopt(a(2), [0 0 0 0 0 b(2)], 1, 1, 1)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% The closed-loop model is calculated using the PIPCL function.
% Finally, both the open and closed loop step responses, together
% with the emissions control input, are plotted using the standard
% Matlab function FILTER.
 
[acl, bcl, bclu]=pipcl(a(2), [0 0 0 0 0 b(2)], v);
 
xd=[zeros(10, 1); ones(90, 1)];    % set-point 
xol=filter([0 0 0 0 0 b], a, xd);  % uncontrolled (open-loop)
x=filter(bcl, acl, xd);            % closed-loop response
u=filter(bclu, acl, xd);           % control input
 
clf; subplot(221)
plot([x xol xd])
subplot(223)
plot(u)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% A faster response can be obtained by increasing the 
% integral-of error weighting to 5000 in the LQ cost function
 
v=pipopt(a(2), [0 0 0 0 0 b(2)], 5000, 1, 1);
 
[acl, bcl, bclu]=pipcl(a(2), [0 0 0 0 0 b(2)], v);
x=filter(bcl, acl, xd);
u=filter(bclu, acl, xd);
 
subplot(222)
plot([x xol xd])
subplot(224)
plot(u)
 
echo off

% end of m-file
