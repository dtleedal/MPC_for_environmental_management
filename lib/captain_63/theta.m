function theta
% THETA  Information about the Captain Toolbox THETA matrix
%
% THETA is a matrix containing information about the transfer
% function model structure, estimated parameters and their
% estimated accuracy. It is generated by a number of Captain
% Toolbox modelling functions. The toolbox function GETPAR is
% generally used to extract the parameter estimates and associated
% covariance matrix.
%
% Example: th=riv([y u], [2 1 3 0]); [a, b]=getpar(th)
%   estimate a transfer function model using the RIV function
%   then extract the denominator and numerator polynomials
%
% See also RIV, RIVID, RIVC, RIVCID, MAR, GETPAR, PREPZ, SCALEB

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

help theta

% end of m-file