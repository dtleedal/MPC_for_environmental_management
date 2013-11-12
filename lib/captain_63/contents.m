% Captain Toolbox
% Version 6.3  23-Dec-2007
%
% Unobserved Components Models.
%   dhr       - Dynamic Harmonic Regression analysis.
%   dhropt    - DHR hyper-parameter estimation.
%   irwsm     - Integrated Random Walk smoothing and decimation.
%   irwsmopt  - IRWSM hyper-parameter estimation.
%   univ      - Trend with Auto-Regression component.
%   univopt   - Trend with AR hyper-parameter estimation.
%
% Time Variable and State Dependent Parameter Models.
%   dar       - Dynamic Auto-Regression and time frequency analysis.
%   daropt    - DAR hyper-parameter estimation.
%   darsp     - DAR spectra plot.
%   darx      - DAR eXogenous variables analysis.
%   darxopt   - DARX hyper-parameter estimation.
%   dlr       - Dynamic Linear Regression analysis.
%   dlropt    - DLR hyper-parameter estimation.
%   dtfm      - Dynamic Transfer Function analysis.
%   dtfmopt   - DTF hyper-parameter estimation.
%   sdp       - State Dependent Parameter analysis.
%
% Transfer Function Models (upgraded: allows for full Box-Jenkins).
%   getparbj - Extract parameters from theta matrix.
%   rivbj    - Discrete-time Transfer Function model estimation.
%   rivbjid  - Discrete-time TF order identification.
%   rivcbj   - Continuous-time TF model estimation.
%   rivcbjid - Continuous-time TF order identification.
%   ivarma   - Noise model estimation.
%   ivarmaid - Noise model identification.
%
% Transfer Function Models (original: for backwards compatibility).
%   getpar    - Extract parameters from theta matrix.
%   riv	      - Discrete-time Transfer Function model estimation.
%   rivid     - Discrete-time TF order identification.
%   rivc      - Continuous-time TF model estimation.
%   rivcid    - Continuous-time TF order identification.
%
% Identification and Diagnostics.
%   acf       - Autocorrelation and Partial Autocorrelation.
%   aic       - Akaike Information Criterium.
%   arspec    - Auto-Regression spectrum.
%   boxcox    - Box-Cox transformation for homoskedasticity.
%   ccf       - Sample Cross-Correlation Function.
%   cusum     - CUSUM and CUSUMSQ tests.
%   histon    - Histogram superimposed over Normal distribution.
%   mar       - Auto-Regresive model estimation.
%   period    - Periodogram estimation.
%   statist   - Sample descriptive statistics.
%
% Univariate control system design.
%   dlqri     - Iterative linear quadratic regulator design.
%   gains     - Proportional-Integral-Plus polynomials.
%   nmssform  - Non-Minimal State Space form.
%   pip       - PIP pole assignment.
%   pipcl     - PIP closed-loop transfer functions.
%   pipcom    - PIP with command input anticipation.
%   piplib    - Simulink library for PIP control.
%   pipopt    - PIP Linear Quadratic optimal.
%
% Multivariate control system design.
%   mfdform   - Matrix Fraction Description form.
%   mfd2nmss  - Multivariable non-minimum state space form.
%   mpipinit  - Initialise Simulink diagram.
%   mpipqr    - Linear Quadratic weightings for PIP control.
%
% Auxiliary Functions.
%   createth  - Creates theta matrix from parameters.
%   del       - Matrix of delayed variables.
%   fcast     - Prepare data for forecasting.
%   kalmanfis - Kalman Filter and Fixed Interval Smoother.
%   mcpar     - Parameters for Monte Carlo analysis.
%   prbs      - Pseudo Random Binary Signal generator.
%   prepz     - Prepare data for input-output modelling.
%   reconst   - Reconstructs a series with jumps.
%   scaleb    - Rescale numerator polynomials.
%   shade     - Plot shaded confidence bounds.
%   stand     - Standardise or de-standardise matrix.
%   theta     - Information about the theta matrix.
%
% Demonstrations.
%   captdemo  - Demonstrations and background information.
%
% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom.
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor.
