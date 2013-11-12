function [F,g,d,h]=mfd2nmss(A,B,inc);
% MFD2NMSS  Multivariable non-minimum state space form
%
% [F,G,D,H]=mfd2nmss(a_mfd,b_mfd,inc)
%
% a_mfd: Matrix Fraction Description (denominator) (*)
% b_mfd: Matrix Fraction Description (numerator) (*)
% inc: Full NMSS (0 - default) or forces the routine to return
%      the regulator NMSS form without integral of error state (1)
%
% F: Transition matrix
% G: Input matrix
% D: Command input matrix
% H: Observation matrix
%
% See also MFDFORM, MPIPQR, MPIPINIT, PIPLIB, DLQRI

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The Matrix Fraction Description, for which amfd represents the
% ‘denominator’ parameters and bmfd the ‘numerator’ parameters, is
% used to determine the NMSS form of a multivariable system.
% Optional inclusion of a positive inc forces the routine to return
% the regulator NMSS form without an integral-of-error state vector.
% The output arguments are the state transition matrix F, input
% matrix G, command input matrix D and observation matrix H.

if nargin==0
  disp(' ')
  disp(' MFD2NMSS  Multivariable non-minimum state space form')
  disp(' ')
  disp(' [F,G,D,H]=mfd2nmss(a_mfd,b_mfd,inc)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, inc=[]; end

[F,g,d,h]=mfd2nmss0(A,B,inc);

% end of m-file