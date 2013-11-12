function [aa,bb]=mcpar(th,nt,ndat,p_level);
% MCPAR  Parameters for Monte Carlo analysis
%
% [aa, bb]=mcpar(th,nt,ndat,p_level)
%
% th: Theta matrix (see 'help theta') (*)
% nt: Number of Monte Carlo realisations (*)
% ndat, plevel: See below.
%
% bb/aa: Each row is one model realisation (trunciated form)
%
% The first row contains non-variated parameters vector
% contained in th. The covariance information from th
% is used in either of two ways:
%
%  1 -  When the routine is called with two input arguments,
%       the pseudo-Gaussian random deviations of parameters
%       are returned
%
%  2 -  When all four arguments are used, the parameters 
%       are generated uniformly over a confidence ellipsoid
%       determined by the covariance matrix eigenvectors and
%       confidence parameter found from the F-distribution
%       with n.d.f. (np, ndat-np) at p-level of confidence
%       If plevel>=1, this parameter is treated as the 
%       number of normalized standard deviations
%
% See also CREATETH

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Here, th is a matrix containing information about the transfer
% function model structure, estimated parameters and their estimated
% accuracy (see theta) and nt is the number of Monte Carlo
% realisations. The function returns, in truncated form, the system
% numerator bb (with assumed unit delay) and denominator aa (with
% assumed leading unity) polynomials, with each row representing one
% Monte Carlo realisation. Typically the function is used to evaluate
% the robustness of Proportional-Integral-Plus (PIP) control systems.
%
% The first row contains the deterministic parameters vector contained
% in th. The covariance information from th is used in either of two
% ways. When the routine is called with two parameters, the pseudo-
% Gaussian random deviations of parameters are returned. Secondly,
% when all four input arguments are used, the parameters are
% generated uniformly over a confidence ellipsoid determined by the
% covariance matrix eigenvectors and confidence parameter found from
% the F-distribution with n.d.f. (nt, ndat-nt) at plevel of confidence.
% If plevel>=1, this parameter is treated as the number of normalized
% standard deviations.

if nargin==0
  disp(' ')
  disp(' MCPAR  Parameters for Monte Carlo analysis')
  disp(' ')
  disp(' [aa,bb]=mcpar(th,nt,ndat,p_level)')
  disp(' ')
  return
end

if nargin<1, th=[]; end
if nargin<2, nt=[]; end
if nargin<3, ndat=[]; end
if nargin<4, p_level=[]; end

if nargin<3
  [aa,bb]=mcpar0(th,nt);
else
  [aa,bb]=mcpar0(th,nt,ndat,p_level);
end

% end of m-file