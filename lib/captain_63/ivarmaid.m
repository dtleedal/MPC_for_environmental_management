 function [p,q,C,D,Pcc,resvar,eef,RR]=ivarmaid(n,nnr,nar,sortcrit,gm)    
% IVARMAID  Identification of ARMA model order for noise signal
%
% [p,q,C,D,Pcc,resvar,eef,RR]=ivarmaid(n,nnr,nar,sortcrit,gm)
%
% n: Noise time series (*)
% nnr: Range of ARMA model orders: e.g. [1 1 4 4] (default [1 1 5 5])
% nar: Order of AR model used to start algorithm (default 50)
% sortcrit: Criterion for sorting: 1=SIC; 2=AIC; 3=FPE
% gm: ARMA estimation algorithm (0)
%        0: IVARMA (CAPTAIN)
%        1: PEM (rather quicker that IVARMA but requires SID Toolbox)
%
% p: AR order 
% q: MA order
% C: AR poly
% D: MA poly
% Pcc: Covariance matrix
% resvar: Residual variance 
% eef: Residuals
% RR=[i,j,AIC,SIC,FPE,cov(eef),R2]
%
% See also IVARMA, RIVBJ, RIVBJID, RIVCBJ, RIVCBJID

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

if nargin==0
  disp(' ')
  disp(' IVARMAID  Identification of ARMA model order for noise signal')
  disp(' ')
  disp(' [p,q,C,D,Pcc,resvar,eef,RR]=ivarmaid(n,nnr,nar,sortcrit,gm)')
  disp(' ')
  return
end

if nargin<1, n=[]; end
if nargin<2, nnr=[]; end
if nargin<3, nar=[]; end
if nargin<4, sortcrit=[]; end
if nargin<5, gm=[]; end

if nargout==1
  p=ivarmaid0(n,nnr,nar,sortcrit,gm);
else
  [p,q,C,D,Pcc,resvar,eef,RR]=ivarmaid0(n,nnr,nar,sortcrit,gm);
end

% end of m-file