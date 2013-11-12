function beraj= histon0(x, nbar, Title)
% HISTON  Histogram over Normal distribution and normality test
%
% beraj=histon(x,nbar,title)
%
% x: Time series (*)
% nbar: Number of bins for histogram (sqrt(length(x)))
% title: Plot title (none)
%
% beraj: Bera-Jarque statistic and P value

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% A histogram of the time series is plotted and the normal 
% theoretical distribution with the same mean and variance is 
% superimposed. This allows for a visual inspection of the 
% distribution of the variable. A formal statistical test, the Bera-
% Jarque test, is also shown, together with its probability value.
% The time series vector y (column) is the only compulsory input 
% to this function. The number of bars to plot in the histogram may 
% be selected by nbar, and Title produces a title in the figure.
% The output is the value of the Bera-Jarque statistic and its 
% probability value.

if nargin==0
  disp(' ')
  disp(' HISTON  Histogram over Normal distribution and normality test')
  disp(' ')
  disp(' beraj=histon(x,nbar,title)')
  disp(' ')
  return
end

if nargin<1, x=[]; end
if nargin<2, nbar=[]; end
if nargin<3, Title=[]; end

beraj=histon0(x, nbar, Title)

% end of m-file