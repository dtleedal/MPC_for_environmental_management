function [x, m, s] = prepz(z, sel, bas, nz, sc, dt)
% PREPZ  Prepare data for input-output modelling
%
% [z,m,s]=prepz(z,sel,bas,nz,sc,dt)
%
% z: I-O data, [Y,U1,...,Unu] where Y and Ui are column vectors (*)
% sel: Select data for modelling [st, en] (1)
%          st: Initial sample
%          en: Final sample (or omit for end)
% bas: Subtract 'baseline' for each data series (1)
%        0: Leave unchanged
%        1: Substract initial value from each series
%        X: Substract mean of first 1:X values
% nz: Prepend each series with 'nz' initial values (0)
% sc: Scale inputs to the same magnitude as y (0 - off)
%       1: Scale each input U1 to Unu
% dt: Subsampling (1)
%
% z: Prepared data in the same format at the input argument
% m: Vector of subtracted baselines
% s: Vector of input gains (see 'help scaleb')
%
% Example: z=prepz([y u1 u2]); th=riv(z, [2 1 1 1 1 0])
%   subtract the initial value from y, u1 and u2, and then
%   estimate a transfer function model using the RIV function
%
% See also RIV, RIVID, RIVC, RIVCID, GETPAR, SCALEB, STAND, THETA

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The only compulsory input argument is z, a matrix of input 
% output data where the first column is considered the output and 
% the rest are inputs.
% 
% The remaining input arguments are optional, with default values 
% set in case they are not specified: sel = [st en] select a sub-
% sample of the original data for modelling from st to en (default 
% is the entire time span); bas tells the function to either subtract 
% the initial value from each series (1 - the default), subtract the 
% mean of the first n values (n), or leave the series unchanged (0); 
% nz adds nz initial values to each series (default is 0); sc sets on 
% or off the scaling of all inputs to the same magnitude as the 
% output (default is 0, off); and dt determines the subsampling 
% interval (default is no subsampling, i.e. 1).
% 
% The output arguments are z, the data prepared according to the 
% input values supplied by the user; m is the vector of subtracted 
% baselines; and s is the vector of input gains. Note that m and s 
% are typically used as input arguments for the function scaleb.

if nargin==0
  disp(' ')
  disp(' PREPZ  Prepare data for input-output modelling')
  disp(' ')
  disp(' [z,m,s]=prepz(z,sel,bas,nz,sc,dt)')
  disp(' ')
  return
end

if nargin<1, z=[]; end
if nargin<2, sel=[]; end
if nargin<3, bas=[]; end
if nargin<4, nz=[]; end
if nargin<5, sc=[]; end
if nargin<6, dt=[]; end

[x, m, s]=prepz0(z, sel, bas, nz, sc, dt);

% end of m-file