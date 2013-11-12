function [CUSUM, CUSUMSQ]= cusum(y)
% CUSUM  CUSUM and CUSUMSQ tests
%
% [CUSUM, CUSUMSQ] = cusum(y)
%
% y: Time series (*)
%
% CUSUM: Cusum recursive test for time varying mean
% CUSUMSQ: Cusumsq recursive test for time varying variance
%
% See also HISTON, BOXCOX

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% This function provides the usual CUSUM (cumulative sum) and 
% CUSUM of squares of a given time series. The CUSUM is 
% computed as the cumulative sum of the standardised series. If 
% the mean changes along the series there will be runs of positive 
% or negative observations; here, the cumulative sum will drift up 
% or down the zero line. The formal assessment is completed on 
% the basis of the 5% confidence bands that are also plotted in the 
% graphical output. Crossing any of the boundaries is an indication 
% of lack of constancy in the mean.
% 
% The CUSUM of squares is similar to the previous description, 
% except that the cumulative sum is performed on the cumulative 
% sum of squares. Therefore, it always increases as time increases, 
% but if the time series is homoskedastic, this increment should be 
% linear, i.e. it should not depart from the 45º line.
% The input y is just the time series, and the outputs CUSUM and 
% CUSUMSQ are the cumulative sum and the cumulative sum of 
% squares, respectively.

if nargin==0
  disp(' ')
  disp(' CUSUM  CUSUM and CUSUMSQ tests')
  disp(' ')
  disp(' [CUSUM, CUSUMSQ] = cusum(y)')
  disp(' ')
  return
end

if nargin<1, y=[]; end

[CUSUM, CUSUMSQ] = cusum0(y);

% end of m-file