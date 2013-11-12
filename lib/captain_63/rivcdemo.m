% RIVCDEMO  Captain Toolbox demonstration
%
% Identification of continuous-time Transfer Function
% models for data from a winding pilot plant
%
% See also RIVC, RIVCID

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% RIVCDEMO  Captain Toolbox demonstration
 
% This script analyses data from a winding pilot plant using the
% continuous time version of the Simplified Refined Instrumental
% Variable (SRIV) algorithm.
 
% The Matlab Control System Toolbox is required for the last
% part of this demo. If this is not installed on your system,
% the script will crash (or press Ctrl+C now to quit).
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The demo uses 3-input, single output data collected from a winding
% process. Winding systems are in general continuous, nonlinear
% processes. They are encountered in a wide variety of industrial
% plants such as rolling mills in the steel industry, plants involving
% web conveyance including coating, papermaking and polymer film
% extrusion processes.
 
% The main role of a winding process is to control the web conveyance
% in order to avoid the effects of friction and sliding, as well as the
% problems of material distortion which can also damage the quality of
% the final product.
 
% The experiments and the data are obtained from:
 
% T. Bastogne, H. Noura, P. Sibille, A. Richard (1998) Multivariable
% identification of winding process by subspace methods for a tension
% control. Control Engineering Practice, vol. 6, no. 9, 1077-1088.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% First load the data.
 
load wind.dat; z=wind;
 
% The system output is:
%   z(:,1) : motor angular speed 2
 
% The system inputs are:
%   z(:,2) : setpoint of motor current 1
%   z(:,3) : setpoint of motor angular speed 2
%   z(:,4) : setpoint of motor current 3
 
Ts=0.01;  % sampling rate (seconds)
N=length(z);  % number of data points
t=[0:0.01:(N-1)*Ts]';  % time vector
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Plot the input data.
 
clf; subplot(3,1,1)
plot(t,z(:,2)); grid; ylabel('Current Setpoint 1');
title('Raw pilot winding input data')
subplot(3,1,2)
plot(t,z(:,3)); grid; ylabel('Speed Setpoint 2');
subplot(3,1,3)
plot(t,z(:,4)); grid; ylabel('Current Setpoint 3');
xlabel('Time in seconds')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause

% Plot the output data.
 
clf; plot(t,z(:,1)); grid;
ylabel('Angular Speed 2');
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Prepare data matrix 'z' from de-meaned data.
 
z=prepz(z, [], length(z));
 
% We will use RIVCID to determine a satisfactory model
% structure, searching for up to 2 parameters each for
% the numerator and denominator polynomials and specifying
% a time delay of zero for each input.
 
nn=[1 1 1 1 0 0 0; 2 2 2 2 0 0 0];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The Simplified Refined Instrumental Variable (SRIV)
% algorithm is employed by using the following settings
% for 'flags' (see 'help rivcid' for details).
 
Ni=8;     % number of SRIV iterations
dt=0.01;  % sampling interval
ddt=1;    % sampling for initial discrete time identification
cf=0;     % adaptive pre-filter flag
Sc=2;     % RT2 selection criterion
flags=[Ni, dt, ddt, cf, Sc];
 
% Finally, we need to specify the prefilter. Here, setting
% the prefilter option to unity instructs the function to estimate
% a discrete time filter and convert this to continous time.
 
c=1;
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% ESTIMATING MODEL : PLEASE WAIT
 
[th,stats,e,rr]=rivcid(z,nn,flags,c);
 
% Scroll up to see the top few models in terms of their RT2
% value, where the highest RT2 (closest to unity) at the top
% of the list indicates the best fit to the data.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Model fit, where unity implies a perfect fit.
 
RT2=stats(3)  % Coefficient of Determination
 
% Extract parameter estimates.
 
[a, b]=getpar(th)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% How good is this model? One way to find out is to plot the model
% output from the returned arguments (measured output-estimated noise)
% and compare this with the measured output.
 
ym=z(:,1)-e;
 
% Set initial model output to correct initial conditions misalignment.
 
ym(1:50)=z(1:50,1);
 
% Let us compare the real and model output.
 
clf
plot(t, [z(:,1) ym z(:,1)-ym+0.2 ones(size(ym))*0.2])
xlabel('Time (seconds)')
title('Data (noisy), simulated and error + 0.2')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We obtain the same plot if we simulate the model using lsim.
 
ys1=lsim(b(1,:), a, z(:,2),t);  % response to input 1
ys2=lsim(b(2,:), a, z(:,3),t);  % response to input 2
ys3=lsim(b(3,:), a, z(:,4),t);  % response to input 3
ys=ys1+ys2+ys3;  % total simulated output
 
ys(1:50)=z(1:50,1);  % initial conditions as before
 
clf
plot(t, [z(:,1) ys z(:,1)-ys+0.2 ones(size(ys))*0.2])
xlabel('Time (seconds)')
title('Data (noisy), simulated and error + 0.2')
 
echo off

% end of m-file