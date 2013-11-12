function y = prbs(N,dt,seed)
% PRBS  Pseudo Random Binary Signal generator
%
% y=prbs(N,dt,J)
%
% N: Required length of the sequence (*)
% dt: Switching period (an integer > 2) (2)
% J: Integer random number generator state (current state)
%
% See also RIV, RIVC

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The PRBS y, is a function of the switching period dt (an integer 
% > 2) and state of the random number generator J (an integer), 
% while N is the length of the series.
    
if nargin==0
  disp(' ')
  disp(' PRBS  Pseudo Random Binary Signal generator')
  disp(' ')
  disp(' y=prbs(N,dt,J)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, dt=[]; end
if nargin<3, seed=[]; end

y=prbs0(N,dt,seed);

% end of m-file