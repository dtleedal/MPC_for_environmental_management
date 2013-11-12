function [inn,zh,xt1,F,Pt1,vr]=kalmanfis(z,u,Phi,E,H,Q,R,Gam,D,P0,x0,sm,Inter,intD)
% KALMANFIS  Fixed Interval Smoother for general state space system
%
% [inn,yhat,xhat,Py,Px,vr]=kalmanfis(y,u,Phi,E,H,Q,R,Gam,D,P0,x0,sm,Int,IntD)
%                                    1 2  3  4 5 6 7  8  9 10 11 12 13   14
% y: Output data (*)
% u: Input data (*)
% E (*), H (*), Q (*), R (1): System matrices (maybe time varying)
% Phi (*), Gam ([]), D ([]): Time invariant system matrices
%   x(t+1)= Phi  x(t) + Gam u(t) + E w(t) ; Q(t) = COV(w(t))
%   y(t)  = H(t) x(t) +  D  u(t) +   v(t) ; R(t) = COV(v(t))
% P0: Initial condition for P matrix (1e5)
% x0: Initial condition for states (0)
% sm: Smoother on (1) or off (0) (default - 1)
% Int: Variance intervention points (none)
% IntD: Intervention type (1e5)
%
% inn:  Innovations
% yhat: Filtered or smoothed fitted values (depending on sm)
% xhat: Estimated states
% Py (F):  Covariance matrix of innovations (filtered or smoothed)
% Px (Pt1): Covariance matrix of states (filtered or smoothed)
% vr: Residual variance
%
% See also DLR, DHR, DAR, DARX, DTFM, IRWSM, SDP

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional Author for KALMANFIS: Renata Romanowicz

% This function gives the user access to quite general Kalman Filter
% (KF) and Fixed Interval Smoothing (FIS) algorithms for the state
% space assimilation and forecasting of uniformly sampled time
% series data. The function is included for an experienced user of
% the toolbox who wishes to access the KF/FIS algorithms directly,
% without resort to the shells for implementing numerous standard
% model types (such as dlr, dhr, dar and darx).
%
% The user provides the data (often input-output data given by u
% and y respectively), as well as the variance/covariance hyper-
% parameters that define the stochastic inputs and observational
% error. Refer to the on-line help for the state space model structure
% and covariance matrices, i.e. Phi, E, H, Q, R, Gam and D. The initial
% state vector and diagonal of the P-matrix may be specified using x0
% and P0, with default values of 0 and 1e5 respectively. FIS may be
% turned off by changing sm from its default unity to 0. In this case,
% the model fit and estimated parameters are their filtered values.
% Int allows for sharp (discontinuous) local changes in the parameters
% at the user supplied intervention points. Here, Int should take the
% same dimensions as y, with positive values indicating variance
% intervention required. Finally, InD gives the diagonal of the
% variance intervention matrix.
%
% The function generates the KF (filtered, forward pass) and FIS
% (smoothed, backward pass) estimates of the state variables xhat and
% output yhat. The innovations series inn, together with the covariance
% matrix of innovations (filtered or smoothed) Py, covariance matrix
% of states (filtered or smoothed) Px and residual variance vr are also
% returned.
%
% Often, the state space model required for this tool will be generated
% by prior RIVID/RIV identification and estimation of a transfer function
% model that is then converted into the required discrete-time state
% space form for kalmanfis. This is illustrated in a demonstration script
% kfisdemo concerned with the assimilation of rainfall-flow data and the
% estimation of unobserved states, where the latter are associated with
% the surface and groundwater flow components that combine to produce the
% measured river flow.

if nargin==0
  disp(' ')
  disp(' KALMANFIS  Fixed Interval Smoother for general state space system')
  disp(' ')
  disp(' [inn,zh,xt1,F,Pt1,vr]=kalmanfis(z,u,Phi,E,H,Q,R,Gam,D,P0,x0,sm,Inter,intD)')
  disp(' ')
  return
end

if nargin<1, z=[]; end
if nargin<2, u=[]; end
if nargin<3, Phi=[]; end
if nargin<4, E=[]; end
if nargin<5, H=[]; end
if nargin<6, Q=[]; end
if nargin<7, R=[]; end
if nargin<8, Gam=[]; end
if nargin<9, D=[]; end
if nargin<10, P0=[]; end
if nargin<11, x0=[]; end
if nargin<12, sm=[]; end
if nargin<13, Inter=[]; end
if nargin<14, intD=[]; end

[inn,zh,xt1,F,Pt1,vr]=kalmanfis0(z,u,Phi,E,H,Q,R,Gam,D,P0,x0,sm,Inter,intD);

% end of m-file