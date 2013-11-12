% RIVDEMO3  Captain Toolbox demonstration
%
% Identification of Transfer Function models from
% simulated input-output data using the graphical interface
%
% See also RIV, RIVID, RIVDEMO1, RIVDEMO2

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% RIVDEMO3  Captain Toolbox demonstration
 
% This script invokes the graphical interface for the
% function RIVID, in order to identify a discrete-time
% transfer function model for a biological data set.
% We start by loading the input-output data.
 
z=load('photo.dat');
plot(z)

% These data show the photosynthetic response of runner
% beans to step changes in ambient carbon dioxide concentration.
% The output variable (first column) is the CO2 assimilation rate
% (umol/m2/s), while the input is the ambient partial pressure of
% CO2 (Pa). These data have already been preprocessed and
% scaled ready for modelling.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% We will use RIVID to determine a satisfactory model
% structure, searching for up to 3 parameters each for
% the numerator and denominator polynomials and for a
% possible 0 to 3 time delays. There is no noise model.
 
nn=[1 1 0 0; 3 3 3 0]
  
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% When the graphical interface is started in a moment,
% A model structure of [1 2 0 0] should appear at the top
% of the list, i.e. first order with 2 numerator parameters,
% no time delay and no noise model. Note that the structure
% with the 2nd lowest YIC value [2 1 0 0] has a relatively
% poor fit, with unrealistic oscillations in the simulated
% response: try clicking in the list box to see this.
 
% However, this preliminary analysis shows that even
% the best fitting linear model does not fully describe
% the nonlinear behaviour of this biological system.
% For this reason, use of DTFM with time variable
% parameters is recommended for the next stage of the
% analysis.
 
% When you are ready to exit the graphical interface,
% simply press the ENTER key (RETURN) and the results will
% be automatically saved to the workspace.
  
% The interface is started by setting the 4th input argument
% in the call to RIVID to unity, as shown below.
 
% --------------------------------------------------------
%          Hit any key to start the user interface
% --------------------------------------------------------
pause

close(1)

th=rivid(z, nn, [], [], 1);  % hit RETURN to continue

% Note that all the standard output arguments from RIVID,
% such as the theta matix (th), are returned to the workspace
% if specified. See 'help rivid' for details. For example:
 
[a, b]=getpar(th)
 
echo off

% end of m-file