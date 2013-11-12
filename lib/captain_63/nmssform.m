function [F,g,d,h]=nmssform(a,b,inc);
% NMSSFORM  Non-minimal state space form
%
% [F,g,d,h]=nmssform(a,b,inc)
%
% b/a: System polynomials (trunciated form) (*)
% inc: If positive function returns regulator form (0)
%
% F: State transition matrix
% g: Input vector
% d: Command vector
% h: Observation matrix

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Truncated form system numerator b (with assumed unit delay) and
% denominator a (with assumed leading unity) polynomials are used
% to determine the NMSS form. Optional inclusion of a positive inc
% forces the routine to return the regulator NMSS form without an
% integral-of-error state vector. The output arguments are the
% state transition matrix F, input vector g, command input vector d
% and observation vector h.

if nargin==0
  disp(' ')
  disp(' NMSSFORM  Non-minimal state space form')
  disp(' ')
  disp(' [F,g,d,h]=nmssform(a,b,inc)')
  disp(' ')
  return
end

if nargin<1, a=[];, end
if nargin<2, b=[]; end
if nargin<3, inc=[]; end

[F,g,d,h]=nmssform0(a,b,inc);

% end of m-file