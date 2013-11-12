% PIPDEMO1  Captain Toolbox demonstration
%
% Univariate Proportional-Integral-Plus (PIP) control
% design for a system with two samples time delay
%
% See also PIPDEMO2, PIPDEMO3

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact

s=warning('off');

echo on

clc
% PIPDEMO1  Captain Toolbox demonstration
 
% This script develops a univariate Proportional-Integral-Plus (PIP)
% control system for a 1st order discrete-time transfer function model
% with two samples pure time delay.
 
at=[1 -0.8];
bt=[0 0 2];
 
% The toolbox control functions assume that all models have at least
% unit delay and include the leading 1 of the demoninator polynomial.
 
a=at(2:end);
b=bt(2:end);

% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% PIP control systems are based on either pole assignment
% or optimisation of a Linear Quadratic (LQ) cost function
% (in the switch below, we have chosen the latter).
 
if 0
  p=[0.6 0.7];                 % two closed-loop poles
  v=pip(a, b, p);              % pole assignment
else
  ew=0.1;                      % control error weight
  uw=1;                        % control input weights
  xw=1;                        % state weights
  v=pipopt(a, b, ew, uw, xw);  % LQ optimal
end
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% PIP controllers include a feedback filter (called fpip)
% including proportional action, an integral controller
% (with gain kpip below) and (assuming the order of the
% numerator polynomial is greater than unity, as here)
% an input filter (called gpip below).
 
[fpip, gpip, kpip]=gains(a, b, v)  % gains
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The closed-loop transfer function and hence step
% command response are found and plotted as follows:
 
[acl, bcl, bclu]=pipcl(a, b, v);
r=zeros(100, 1);         % command input
r(10:59)=ones(50, 1);    % insert step
y=filter(bcl, acl, r);   % output variable
u=filter(bclu, acl, r);  % input variable
 
subplot(211)
plot([y r])
axis([0 length(r) -0.2 1.2])
title('Closed-loop output')
subplot(212)
plot(u)
axis([0 length(r) -0.05 0.14])
title('Input signal')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Monte Carlo Simulation (MCS) is used to investigate the
% robustness of the above controller to model mismatch
 
% In the first case, RIV (see 'help riv') is utilised to
% generate parametric uncertainty with a realistic
% covariance distribution for this model structure
 
u=[zeros(10, 1); prbs(90, 5, 0)];    % input signal
y=filter(bt, at, u)+randn(size(u));  % simulated output
th=riv([y u], [1 1 2 0]);            % estimate model
[atm, btm, tmp, P]=getpar(th);       % extract parameters
ym=filter(btm, atm, u);              % model response
clf; plot([y ym u])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Create matrix of parameters for 50 Monte Carlo realisations
[aa, bb]=mcpar(th, 50);
 
% The loop below determines the closed-loop transfer function
% for each realisation and plots the associated response and
% pole positions. Two PIP control structures are compared,
% namely the standard feedback form (FB) and an alternative
% forward path (FP) form. The latter structure uses an internal model,
% hence the use of 'a' and 'b' as input arguments in the call
% to PIPCL below.
 
% for ff=1:length(aa)
%   [acl, bcl]=pipcl(aa(ff, :), bb(ff, :), v);  % FB form
%   rfb(:, ff)=roots(acl);                      % poles
%   yfb(:, ff)=filter(bcl, acl, r);             % response
%   [acl, bcl]=pipcl(aa(ff, :), bb(ff, :), v, a, b);  % FP form
%   rfp(:, ff)=roots(acl);                            % poles
%   yfp(:, ff)=filter(bcl, acl, r);                   % response
% end
 
% --------------------------------------------------------
%  Hit any key to activate this loop and plot the results
% --------------------------------------------------------
pause

echo off
for ff=1:length(aa)
  [acl, bcl]=pipcl(aa(ff, :), bb(ff, :), v);
  rfb(:, ff)=roots(acl);
  yfb(:, ff)=filter(bcl, acl, r);
  [acl, bcl]=pipcl(aa(ff, :), bb(ff, :), v, a, b);
  rfp(:, ff)=roots(acl);
  yfp(:, ff)=filter(bcl, acl, r);
end
echo on

clf
subplot(221)
plot(yfb, 'b')
axis([0 length(r) -0.2 1.2])
title('Feedback response')
subplot(222)
plot(real(rfb), imag(rfb), 'b.')
title('Feedback poles')
axis([-1 1 -1 1])
subplot(223)
plot(yfp, 'b')
axis([0 length(r) -0.2 1.2])
title('Forward path response')
subplot(224)
plot(real(rfp), imag(rfp), 'b.')
title('Forward path poles')
axis([-1 1 -1 1])

% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% If the Simulink extension to Matlab is installed on your
% system, then it can be convenient to use this package
% rather than PIPCL and FILTER for closed-loop simulation.
 
% Once this script has finished, the Simulink diagram can
% be opened by typing 'delaypip' at the command line.
 
% Note that if Simulink is not installed, the last part of
% this script will crash (or press Ctrl+C now to quit).
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% In the loop below, the Simulink diagram 'delaypip' is called
% for each Monte Carlo realisation. In this example, a load
% disturbance is introduced at the 30th sample.
 
% for ff=1:50
%   as=[1 aa(ff, :)];
%   bs=[0 bb(ff, :)];
%   sim('delaypip')
%   yyfb(:, ff)=yfb;
%   yyfp(:, ff)=yfp;
%   uufb(:, ff)=ufb;
%   uufp(:, ff)=ufp;
% end
 
% --------------------------------------------------------
%  Hit any key to activate this loop and plot the results
% --------------------------------------------------------
pause
 
t=[1:length(r)]';  % time vector
 
% ** Please wait for Monte Carlo Simulation to finish **
 
echo off
for ff=1:length(aa)
  as=[1 aa(ff, :)];
  bs=[0 bb(ff, :)];
  sim('delaypip')
  yyfb(:, ff)=yfb;
  yyfp(:, ff)=yfp;
  uufb(:, ff)=ufb;
  uufp(:, ff)=ufp;
end
echo on

clf
subplot(211)
plot(yyfb, 'b')
hold on
plot(r)
axis([0 length(r) -0.2 1.4])
title('Feedback response')
subplot(212)
plot(yyfp, 'b')
hold on
plot(r)
axis([0 length(r) -0.2 1.4])
title('Forward path response')

echo off

warning(s);

% end of m-file