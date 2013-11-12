function [th,e,y0,AH,PH] = mar(y,na,a0,P0)
% MAR  Univariate Auto-Regressive model estimation using Least Squares
%
% [th,e,y0]=mar(y,na,a0,P0)
%
% y: Univariate data series (*)
% na: AR order definition (*)
%       Scalar: Full model
%       Vector: Subset AR
%         e.g. na=[1 2 11:14] implies a3,...,a10 = 0
% a0: Initial parameter vector, presence of this argument
%     forces recursive calculation, non-recursive by default ([])
% P0: initial covariance matrix for recursive calculations (100*I)
%
% th: Theta matrix (see 'help theta')
% e: Error series; model values are equal to y-e
% y0: Interpolated data
% AH: Recursive parameter estimates
% PH: Recursive P-matrix estimates
%
% See also AIC, GETPAR, THETA, UNIV, UNIVOPT

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user. The 
% function automatically handles missing values in y by switching 
% to recursive mode. In fact, y may be appended with additional 
% NaNs to forecast or backcast beyond the original series. The AR 
% model structure is defined by na, which is a scalar or vector 
% listing the required past output variables used in the model. For 
% example, [1:5, 20] specifies a model based on y(t-1) to y(t-5) plus 
% a y(t-20) component (i.e. subset AR). a0 and P0 are the optional 
% initial conditions for the parameters and covariance matrix 
% respectively, a vector of zeros and an identity matrix with 100 
% diagonals by default. Note that specifying a 3rd input argument 
% forces a recursive calculation, rather than the default en block 
% solution.
% 
% The output argument th contains information about the model 
% structure, estimated parameters and their estimated accuracy. 
% The toolbox function getpar is generally used to extract the 
% parameter estimates and associated covariance matrix. The error 
% series e are defined by the model response subtracted from the 
% data, i.e. the model response may be found from y - e. Finally, 
% the interpolated data y0 consist of the original series with any 
% missing data replaced by the model.

if nargin==0
  disp(' ')
  disp(' Univariate Auto-Regressive model estimation using Least Squares')
  disp(' ')
  disp(' [th,e,y0]=mar(y,na,a0,P0)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, na=[]; end
if nargin<3, a0=[]; end
if nargin<4, P0=[]; end

[th,e,y0,AH,PH]=mar0(y,na,a0,P0);

% end of m-file