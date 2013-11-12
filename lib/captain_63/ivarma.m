function [C,D,Pcc,resvar,eef,R2,AIC,SIC,FPE,AR,AH,Pse1,PH]=...
                                          ivarma(n,nn,arorder,nit,rec)
% IVARMA  Estimation of ARMA models for noise signal
%
% [C,D,Pn,resvar,eef,R2,AIC,SIC,FPE,AR,AH,AHse,PH]=ivarma(n,nn,nar,nit,rec)
%
% Or for SID object code output (SID Toolbox required)
%   M = ivarma(n,nn,arorder,nit,rec)
%
% n: Noise time series (*)
% nn: ARMA model order: e.g. [2 1] (*)
% nar: Order of AR model used to start algorithm (default AIC or 50)
% nit: Number of iterations (6)
% rec: Recursive (1) or en bloc (0)
%
% C: ARMA denominator polynomial when called with multiple output arguments
%    If only one output argument: provides model structure of SID type
%    (if SID Toolbox exists) else a structure with the same fields as SID
% D: ARMA numerator polynomial
% Pn: Covariance matrix
% resvar: Residiual variance
% eef: Residual series
% R2: Coefficient of Determination
% AIC: Akaike Inf. Criterion
% SIC: Swartz Inf. Criterion
% FPE: Final Prediction Error
% AR: High oder AR model
% AH: Recursive parameter estimates
% AHse: Standard errors on AH
% PH: Covariance Matrices for AH
%
% See also IVARMAID, RIVBJ, RIVBJID, RIVCBJ, RIVCBJID

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

if nargin==0
  disp(' ')
  disp(' IVARMA  Estimation of ARMA models for noise signal')
  disp(' ')
  disp(' [C,D,Pn,resvar,eef,R2,AIC,SIC,FPE,AR,AH,AHse,PH]=ivarma(n,nn,nar,nit,rec)')
  disp(' ')
  return
end

if nargin<1, n=[]; end
if nargin<2, nn=[]; end
if nargin<3, arorder=[]; end
if nargin<4, nit=[]; end
if nargin<5, rec=[]; end

if nargout==1
  C=ivarma0(n,nn,arorder,nit,rec);
else
  [C,D,Pcc,resvar,eef,R2,AIC,SIC,FPE,AR,AH,Pse1,PH]=ivarma0(n,nn,arorder,nit,rec);
end

% end of m-file