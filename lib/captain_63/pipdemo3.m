% PIPDEMO3  Captain Toolbox demonstration
%
% Multivariate Proportional-Integral-Plus (PIP) control design
%
% See also PIPDEMO1, PIPDEMO2

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% PIPDEMO3  Captain Toolbox demonstration
 
% This script develops a multivariate Proportional-Integral-Plus (PIP)
% control system for a discrete-time model of a coupled drives system.
 
% The Simulink extension to Matlab and the Control System Toolbox
% are required for this demo. If these are not installed on your
% system, the script will crash (or press Ctrl+C now to quit).
 
% Note that the toolbox control functions can be used without
% these extensions to Matlab. However, they are useful for plotting
% and simulating the results, hence are used in this demo.
 
% The system has two inputs: the voltage applied to two motors.
% It has two outputs: the tension and speed of an elastic band.
 
% --------------------------------------------------------
%             Hit any key to open Simulink diagram
% --------------------------------------------------------
pause
 
% The model DRIVEOPEN has been created to determine the
% open loop response of the system. Open the "coupled drives"
% subgroup to see the 4 transfer functions, representing each
% of the input output pathways. The open loop response is shown
% in the scope windows: one for the inputs, one for the outputs.
 
driveopen; sim('driveopen');
 
% In response to a step in each motor voltage, the speed (yellow)
% follows a first order response, while the tension (magenta)
% has an oscillatory second order response. Also, the gain of
% the tension model is positive for motor 1 and negative for
% motor 2. Finally, it is clear that these outputs are highly
% coupled, with both outputs responding to a step in each motor.
 
% ** Please CLOSE this Simulink model before continuing. **
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% For the purposes of this demo, we will assume that the
% system model is unknown. RIVID is used to identify multiple
% input, single output transfer functions for the speed and
% tension, using the data collected from the Simulink scopes.
 
% See 'help rivid' or 'rivdemo' for more information about
% system identification using the Captain Toolbox.
 
% Speed model
z=[yy(:, 2) uu(:, 2:3)];
nn=[1 1 1 1 1 0; 2 2 2 1 1 0];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
[th, stats, e] = rivid(z, nn);
[at_speed, bt_speed]=getpar(th);
subplot(211)
plot(z(:, 1), '.')
hold on
plot(z(:, 1)-e)
axis([0 length(z) -0.1 0.6])
title('Speed')
 
% Tension model
z=[yy(:, 3) uu(:, 2:3)];
nn=[1 1 1 1 1 0; 2 2 2 1 1 0];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
[th,stats,e] = rivid(z,nn);
[at_tension, bt_tension]=getpar(th)
subplot(212)
plot(z(:, 1), '.')
hold on
plot(z(:, 1)-e)
axis([0 length(z) -0.6 0.6])
title('Tension')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% As would be expected for a deterministic simulation, the
% models identified have the same parameters as those in
% the original Simulink model. Note the common denominators
% for the speed (at_speed) and tension (at_tension). By
% contrast, the numerators for each model are represented
% in matrix form, with the first row representing the response
% to input 1 and the second row input 2. The overall response
% would be found by adding these two components.
 
disp(at_speed)
disp(bt_speed)
disp(at_tension)
disp(bt_tension)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% It is useful to save these parameters in two matrices,
% one for the numerator (b_all), the other for the denominator
% (a_all). Here, each row of the matrix represents one of the
% input (u) output (y) pathways, listed in the following order:
% u(1)->y(1), u(2)->y(1), u(1)->y(2), u(2)->y(2). Finally, we
% remove the leading unity and assumed unit delay.
 
a_all=[at_speed 0; at_tension];
a_all=kron(a_all, ones(2, 1));
b_all=[bt_speed zeros(2, 1); bt_tension];
a_all=a_all(:, 2:end)
b_all=b_all(:, 2:end)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Moving onto control system design, we determine the left
% matrix fraction description and associated non-minimal
% state space form as follows:
 
[A_mfd, B_mfd]=mfdform(a_all, b_all);
[F_nmss, G_nmss, D_nmss, H_nmss]=mfd2nmss(A_mfd, B_mfd, 0);
 
% In this example, we assign unity weights to the sum of the
% output states, the sum of the input states and to each
% integral of error state variable, before computing the
% control gains that minimise the standard Linear Quadratic
% cost function for multivariable PIP control.
 
yw=[1 1];  % output weights
uw=[1 1];  % input weights
zw=[1 1];  % IOE weights
[Q, R]=mpipqr(A_mfd, B_mfd, zw, uw, yw);
v=dlqri(F_nmss, G_nmss, Q, R);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The theoretical step response is calculated using the
% state space model in closed-loop form.
 
F_closed=(F_nmss-G_nmss*v);
D_closed=zeros(size(H_nmss, 1), size(D_nmss, 2));
sys=ss(F_closed, D_nmss, H_nmss, D_closed, -1);
y=step(sys, [1:1:160]);
y=[y(:,1,1) y(:,2,1) y(:,1,2) y(:,2,2)];

% for i = 1:4
%   eval(['subplot(2,2,' int2str(i) ');'])
%   plot(y(:,i)); set(gca, 'ylim', [-0.1 1.2])
% end
 
% --------------------------------------------------------
%  Hit any key to activate this loop and plot the results
% --------------------------------------------------------
pause

echo off
clf
for i = 1:4
  eval(['subplot(2,2,' int2str(i) ');'])
  plot(y(:,i)); set(gca, 'ylim', [-0.1 1.2])
  set(gca, 'xlim', [0 length(y)])
end
echo on
 
% The off-diagonal subplots show the cross coupling
% terms, whilst the diagonals show a closed-loop step
% in the speed (top left) and tension (bottom right)
 
% The toolbox function MPIPINIT initialises the Simulink
% block library for multivariable PIP control. The first
% output argument is required for the feedback control
% structure, the second for the forward path structure.
 
[pipfb, pipfp]=mpipinit(A_mfd, B_mfd, v);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% response to a step in speed using Simulink
level=[1 0]; sim('drivepip');

clf; subplot(211)
plot(yy(:, 1), yy(:, [2 4]), 'b', yy(:, 1), yy(:, [3 5]), 'r')
set(gca, 'ylim', [-0.2 1.2])
title('Outputs: Speed (b), Tension (r)')
subplot(212)
plot(yy(:, 1), uu(:, 2), 'b', yy(:, 1), uu(:, 3), 'r')
title('Inputs: Motor 1 (b), Motor 2 (r)')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% response to a step in tension using Simulink
level=[0 1]; sim('drivepip');
 
clf; subplot(211)
plot(yy(:, 1), yy(:, [2 4]), 'b', yy(:, 1), yy(:, [3 5]), 'r')
set(gca, 'ylim', [-0.2 1.2])
title('Outputs: Speed (b), Tension (r)')
subplot(212)
plot(yy(:, 1), uu(:, 2), 'b', yy(:, 1), uu(:, 3), 'r')
title('Inputs: Motor 1 (b), Motor 2 (r)')
 
% The Simulink diagram can now be examined by typing
% 'drivepip' at the command line. This simulation uses the
% feedback form of PIP control. The alternative forward
% path structure can be used by copying the PIP-MFP block
% from the Simulink library 'piplib'.
 
echo off

% end of m-file