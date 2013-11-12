% UNIVDEMO  Captain Toolbox demonstration
%
% Trend and Auto-Regression modelling
% of the air passenger time series
%
% See also UNIV, UNIVOPT, DHRDEMO

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact
echo on

clc
% UNIVDEMO  Captain Toolbox demonstration
  
% This script analyses the air passenger time series using a
% model based on a trend plus an Auto-Regression (AR) component.
 
load air.dat;  % thousands of passengers per month (1949-1960)
 
% Missing values are added in order to test the ability of the
% algorithms to automatically handle interpolation and
% forecasting. Here, fcast is employed to replace samples
% 84-94 (interpolation) and 131-144 (forecasting) with
% Not-a-Number variables.
 
y=fcast(air, [84 94; 131 144]);
 
clf; plot([air y])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% In the first place, we find a smooth trend using an
% Integrated Random Walk (IRW) with an initial Noise
% Variance Ratio (NVR) hyper-parameter of 1e-4.
 
TVP=1;
nvr0=1e-4;
[fit, fitse, trend]=univ(y, [], TVP, nvr0);
clf; plot([y trend y-trend])
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The order of the AR model for the perturbations
% is identified using the Akaike Information
% Criterium.
 
% We will search for models ranging from 1st order
% through to 20th order.
 
% CALCULATING AIC : PLEASE WAIT
 
p=aic(y-trend, 20);  % p is the best fit AR polynomial
 
p=length(p)-1  % order of the AR model
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The AR model for the peturbations is estimated, together
% with the associated optimal NVR for the IRW trend.
 
[nvr, ARp]= univopt(y, [1:p], TVP, nvr0);
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% Finally, both components are estimated simultaneously.
 
[fit,fitse,trend,trendse,comp,y0]=univ(y, ARp, TVP, nvr);
 
clf
plot([trend air])
title('Data and trend')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The seasonal component increases over time.
 
clf
plot(comp)
set(gca, 'xlim', [0 length(air)])
title('Total seasonal component')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% The functions successfully interpolate over the 10 months
% of missing data starting at sample 85. Similarly, they predict
% the output for the final year, i.e. sample 132 until the end
% of the series. Note that since missing data were introduced at
% the start of the analysis, the latter is a 'true' forecast.
 
clf
plot([fit air])
hold on
plot([fit-air zeros(144, 1)])
plot([84, 84], [-50 50])  % start of interpolation
plot([94, 94], [-50 50])  % end of interpolation
plot([131, 131], [-50 50])  % forecasting horizon
set(gca, 'xlim', [0 length(air)])
title('Data, fit and residuals')
 
% --------------------------------------------------------
%                 Hit any key to continue
% --------------------------------------------------------
pause
 
% A zoomed in view of the final two years are shown, with
% the standard errors and forecasting horizon also indicated.
 
clf
plot(fit)
hold on
plot(air, 'o')
plot([fit+2*fitse fit-2*fitse], ':r')
plot([131, 131], [200 700])  % forecasting horizon
axis([119 length(fit) 250 700])
title('Data (o), fit and standard errors')
 
% Note that Integrated and Double Integrated AR models
% for the trend are also possible using univ and univopt,
% so that the full model may take one of the following
% forms: RW+AR; IRW+AR (as above); IAR+AR; and DIAR+AR.
 
echo off

% end of m-file