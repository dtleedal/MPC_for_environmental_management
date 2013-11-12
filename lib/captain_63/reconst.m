function yr=reconst0(y, Interv)
% RECONST  Reconstructs a series with jumps at intervention points
%
% yr=reconst(y,Int)
%
% y: Time series (*)
% Int: Vector of variance intervention or jump points (*)
%
% yr: Reconstructed series without jumps
%
% See also IRWSM, DHR

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The purpose of this function is to reconstruct the behaviour of a 
% series under the hypothesis that several given sudden jumps 
% would never have happened.
% 
% The first input argument y is the original time series, while Int is 
% a vector of intervention points at which the jumps are observed. 
% The output argument yr is the modified time series.

if nargin==0
  disp(' ')
  disp(' RECONST  Reconstructs a series with jumps at intervention points')
  disp(' ')
  disp(' yr=reconst(y,Int)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, Interv=[]; end

yr=reconst0(y, Interv);

% end of m-file