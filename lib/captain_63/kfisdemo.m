% KFISDEMO  Captain Toolbox demonstration
%
% Rainfall-flow analysis and forecasting example
% based on data from the ephemeral Canning River
% in Western Australia
%
% See also KALMANFIS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for KFISDEMO: Renata Romanowicz

echo off
clear all
close all
format compact

y=load('canningflow.dat'); yr=y;
u=load('canningrain.dat'); ur=u;
T=load('canningtemp.dat'); evap=T;

% data start 01/01/1977
% analysis data set
y1=yr(3000:3700);
u1=ur(3000:3700);
evapr=T(3000:3700);
evapmr=evapr-mean(evapr);
evapmr=evapmr';
t=load('canningtime.dat')';
tt=t(3000:3700);
nt=length(y1);
ax=[min(tt), max(tt), 0, 4.5];

clc
echo on
% KFISDEMO  Captain Toolbox demonstration
 
% This script uses the function KALMANFIS for
% rainfall-flow analysis and forecasting
% based on data from the ephemeral Canning River
% in Western Australia
 
% Note that this demonstration is designed for an
% experienced user of the toolbox who is already
% familiar with the methods involved. Please see,
% e.g. DLRDEMO for a more straightforward example,
% using a readily accessible shell for the call to
% Kalman Filtering and Fixed Interval Smoothing.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Plot of data: flow, rainfall and temperature
 
subplot(311)
plot(t(3000:3700), y1, 'k')
axis([1985.235, 1987.14, 0, 4])
title('Flow, Rainfall and Temperature: Canning River, W.A., 1985.2-1987.1')
ylabel('Flow (cumecs)')
subplot(312)
plot(t(3000:3700), u1, 'k')
axis([1985.235, 1987.14, 0, 80])
ylabel('Rainfall (mm)')
subplot(313)
plot(t(3000:3700),evapr,'k')
axis([1985.235, 1987.14, 0, 40])
xlabel('Date')
ylabel('Temperature (deg.C)')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Although not shown here, the RIVID estimation tool was used to
% identify an initial linear Transfer Function (TF) model between
% the effective rainfall and the flow data: here, the effective
% rainfall is obtained as a power law in the flow variable (y1),
% where the flow variable is acting as a surrogate measure of the
% catchment water storage (soil moisture) and the power law
% exponent (pval) has been obtained by prior State Dependent
% Parameter (SDP) model identification and estimation. The
% identified linear TF model has structure [2 2 1 0] where,
% initially, a noise model has not been included for simplicity.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% parameters obtained from SDP and TF analysis
 
pval=0.7776;
mo=[2 2 1 0];
uf=u1.*(y1.^real(pval));
uf=uf.*(sum(y1)/sum(uf));
yf=y1;
Z1=[yf uf];
[TH1, stats, e, var, Ps, Pc, y0]=riv(Z1, mo, [4 2 0 0 0]);
[A1, B1, C, P]=getpar(TH1);
ym1=mpefilt(B1, A1, uf, yf(1:10));
[R, P, K] = residue(B1, A1);
A11=[1 -P(1)];
B11=R(1);
A22=[1 -P(2)];
B22=R(2);
F=eye(2); F(1, 1)=P(1); F(2, 2)=P(2);
B(1, 1)=R(1); B(2, 1)=R(2);
G=[1 1]; D=0; He=[1 1];
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% KALMAN FILTERING
 
% The optimised Noise Variance Ratio (NVR) parameters (Q and R1)
% required to run KALMANFIS have been obtained in prior analysis
% for the above TF model (without a noise model). KALMANFIS is
% now run to generate the flow estimate and the estimates of the
% two state variables, which represent the quick (surface processes)
% and slow (groundwater processes) components of flow. Note: these
% are filtered estimates, hence sm=0 in the function call.
 
Q=[1.0e-003*0.1959 0; 0 1.0e-003*0.0231];  % Q optimised / guessed
sm=0;
R1=0.001*cov(y1)*ones(length(y1),1);
[inn1,yhat,xhat,Py,Px,vr]=kalmanfis(y1,uf,F,G,He,Q,R1,B,D,[],[],sm);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% model fit
rt1=1-cov(inn1)/cov(y1)

% plot results
clf
zf=shade(tt', yhat, sqrt(vr),y1);
axis(ax)
title('DBM model estimation: 1985.2-1987.1: 96.2% of flow explained')
ylabel('Flow (cumecs)')
xlabel('Date')

% This plot shows the KALMANFIS generated Kalman filter prediction of
% the flow, with 96.2% of the flow explained by the predictions.
% Here and subsequently, the shaded area is the 95% confidence
% interval (2 x standard deviation) on the predictions.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
nx=mo(1);
ax=axis;
echo off
px=reshape(Px', nx, nx, nt);
for i=1:nt
    pp=squeeze(px(:,:,i));
    px1(i)=pp(1, 1);
    px2(i)=pp(2, 2);
end
echo on
 
clf
subplot(211)
zf=shade(tt', xhat(:, 1), sqrt(px1'), xhat(:, 1));
axis([min(tt), max(tt), 0, 4])
xlabel('Date')
ylabel('Slow component of runoff')
title('Unobserved components without smoothing')
subplot(212)
zf=shade(tt', xhat(:, 2), sqrt(px2'), xhat(:, 2));
axis([min(tt), max(tt), 0, 4])
xlabel('Date')
ylabel('Fast component of runoff')

% These plots show the KALMANFIS generated Kalman filter 
% estimates of the slow and fast components of the flow.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% FIXED INTERVAL SMOOTHING
 
% Now KALMANFIS is run again with sm=1, to provide the smoothed
% estimates of the flow and the two flow components.
  
sm=1;
[inn2,yhat,xhat,Py,Px,vr]=kalmanfis(y1,uf,F,G,He,Q,R1,B,D,[],[],sm);

% model fit
rt2=1-cov(inn2)/cov(y1)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
clf
zf=shade(tt', yhat, sqrt(vr), y1);
axis(ax)
title('DBM model with smoothing: 1985.2-1987.1: 99.84% of flow explained')
ylabel('Flow (cumecs)')
xlabel('Date')
 
% This plot shows the KALMANFIS generated FIS estimates of
% the flow, with 99.84% of the flow explained by the 
% estimated output.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
nx=mo(1);
px=reshape(Px', nx, nx, nt);
echo off
for i=1:nt
    pp=squeeze(px(:, :, i));
    px1(i)=pp(1, 1);
    px2(i)=pp(2, 2);
end
echo on
 
clf
subplot(211)
zf=shade(tt', xhat(:, 1), sqrt(px1'), xhat(:, 1));
axis([min(tt), max(tt), 0,4])
xlabel('Date')
ylabel('Slow component of runoff')
title('Unobserved components after smoothing')
subplot(212)
zf=shade(tt', xhat(:, 2), sqrt(px2'), xhat(:, 2));
axis([min(tt), max(tt), 0, 4])
xlabel('Date')
ylabel('Fast component of runoff')

% These plots show the KALMANFIS generated FIS estimates 
% of the slow and fast components of the flow. The next step
% is to identify an AR model for the residuals using AIC
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
%
pause

maic=aic(inn1, 30, 1);
 
% estimate AR(25) noise model
TH=mar(inn1, 25);
[a1, Bn, Cn, P, d]=getpar(TH);
 
% estimate white noise input
erm=filter(a1, 1, inn1);
 
% RT2 of the noise model
1-cov(erm)/cov(inn1)
 
% Here, the residuals are modelled as an AR(25) noise process
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The new state space model, incorporating the noise model,
% is now defined - open the script KFISDEMO to see the commands
% for generating this model
echo off

ne=length(a1)-1;
m=length(F);
nc=m+ne;
nu=1;
B1=zeros(nc, nu);
A1=zeros(nc);
% He=zeros(nc,1)';
C1=zeros(nc, 1)';  % instead of He
C1(1:m)=He;
C1(end)=1;
ann=fliplr(a1);
for i=1:m
    for j=1:m
        A1(i, j)=F(i, j);
    end
end
for i=1:ne
    A1(end, m+i)=-ann(i);
end
for i=1:ne-1
    A1(i+m, i+m+1)=1;
end
for i=1:m
    B1(i,1)=B(i);
end
% B(mo(1)+1)=1;
% only mo(1) state vectors have noise + one noise state: mo+1
G1=zeros(nc, m+1);
for i=1:m
    G1(i, i)=1;
end
G1(end, end)=1;
% G=[1 0 0; 0 1 0; 0 0 0; 0 0 0 ; 0 0 0; 0 0 1];  % for the noise
% P=F*P*F' + G*Q*G';
% state noise is only in mo(1)+1 states
% G: [nc,mo(1)+1]
% xk+1=F*xk+B*u+G*Q
D=0;
% y=He*xk+D*uk+R
% He=[1 1 0 0 0 1]
Q1=zeros(m+1);
for i=1:m
    for j=1:m
        Q1(i, j)=Q(i, j);
    end
end
Q1(m+1, m+1)=0.1*cov(erm);  % first estimate of cov(er)
nt=length(y1);
R1=cov(y1)*ones(length(y1), 1)*0.001;

% In practice, the NVR hyper-parameters are found from a numerical
% optimisation step as follows:
%   x0=[-8.8477 -6.6725 -7.0051];
%   [x,fval,exitflag,output]=...
%     fminsearch('kalm_opt52',x0,options,y,uf,uf,A1,G1,C1,R1,B1,D,mo,...
%     np,padapt,hetero,[],sm,noise);
% where kalmopt52 is the user defined goal function (not shown here)
echo on
 
% For this demonstration, we will use the following previously
% optimised values for the NVR hyper-parameters
x=[-8.8477 -6.6725 -7.0051];
 
Q1=zeros(m+1);
for i=1:m+1
  Q1(i, i)=exp(x(i));
end 
  
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
%
pause
 
% The new state space model is now used for Kalman filtering
sm=0;
[inn,yhat,xhat,Py,Px,vr]= kalmanfis(y1,uf,A1,G1,C1,Q1,R1,B1,D,[],[],sm);

% model fit
1-cov(inn)/cov(y1)
 
clf
zf=shade(tt', yhat, sqrt(vr), y1);
axis(ax)
title('DBM + noise Model: 1985.2-1987.1: 96.88% of flow explained')
ylabel('Flow (cumecs)');
xlabel('Date')
 
% This plot shows the KALMANFIS generated Kalman filtering 
% results after the optimisation of the NVR parameters for 
% the state space model incorporating the noise model, with
% 96.9% of the flow explained by the Kalman Filter predictions.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
%
pause
 
maic=aic(inn, 30, 1);
 
clf
acf(inn, 20, [], 1);
 
% This plot shows the autocorrelation function of the final
% residuals (one-step-ahead prediction errors) based on the
% AR(25) noise model.
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
%
pause
 
% VALIDATION
%
% An important part of rainfall-flow modelling is
% predictive validation of the model on data over
% other years that have not been used in the model
% estimation, as shown below.
 
yt=yr(1:500);
ut=ur(1:500);
tt=t(1:500);
uft=ut.*(yt.^real(pval));
uft=uft.*(sum(yt)/sum(uft));
R1=0.001*cov(y1)*ones(length(yt),1);
sm=0;
[inn,yhat,xhat,Py,Px,vr] = kalmanfis(yt,uft,A1,G1,C1,Q1,R1,B1,D,[],[],sm);
Rt2t=1-cov(inn)/cov(yt)
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
%
pause
 
clf
zf=shade(tt',yhat,sqrt(vr),yt);
title('Predictive validation: 1977-1978.5: 93.59% of flow explained')
xlabel('Date')
ylabel('Flow (cumecs)')
axis([1977.4,1978.1,0,3.5]);
 
% This plot shows the predictive validation results produced 
% by KALMANFIS, where 93.60% of the data is explained by the 
% KALMANFIS generated Kalman Filter predictions.
 
echo off

% end of m-file
