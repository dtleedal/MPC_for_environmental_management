function [pipfb, pipfp]=mpipinit(A_mfd, B_mfd, v)
% MPIPINIT  Initialise block diagram for multivarible PIP control
%
% [pipfb,pipfp]=mpipinit(a_mfd,b_mfd,k)
%
% a_mfd: Matrix Fraction Description (denominator) (*)
% b_mfd: Matrix Fraction Description (numerator) (*)
% k: Control gain matrix (*)
%
% pipfb: Cell with feedback PIP settings (for PIPLIB.MDL)
% pipfp: Cell with forward path PIP settings (for PIPLIB.MDL)
%
% See also MFD2NMSS, MPIPQR, MPIPQR, PIPLIB, DLQRI

% Returns MATLAB cells for initialising both feedback pipbf and
% forward path pipfp multivariable PIP control blocks in piplib
% (SIMULINK library). The input arguments are the Matrix Fraction
% Description, for which amfd represents the 'denominator'
% parameters and bmfd the 'numerator' parameters, together with
% the control gain matrix k.

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

if nargin==0
  disp(' ')
  disp(' MPIPINIT  Initialise block diagram for multivarible PIP control')
  disp(' ')
  disp(' [pipfb,pipfp]=mpipinit(a_mfd,b_mfd,v)')
  disp(' ')
  return
end

if nargin<1, A_mfd=[]; end
if nargin<2, B_mfd=[]; end
if nargin<3, v=[]; end

[pipfb, pipfp]=mpipinit0(A_mfd, B_mfd, v);

% end of m-file