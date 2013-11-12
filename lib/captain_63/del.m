function yd=del(y,n)
% DEL  Matrix of delayed variables
%
% yd=del(y,n)
%
% y: time series (*)
% n: maximum delay (1)
%
% yd: nxm matrix of delayed variables for i=1:n

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% The time series vector y (column) is specified by the user, while 
% the optional parameter n (default unity) is the maximum lag 
% required. The output yd is a matrix with n columns, where the 
% first column is y lagged by one sample, the second column is y 
% lagged by two samples and so on.

if nargin==0
  disp(' ')
  disp(' DEL  Matrix of delayed variables')
  disp(' ')
  disp(' yd=del(y,n)')
  disp(' ')
  return
end

if nargin<1, y=[]; end
if nargin<2, n=[]; end

yd=del0(y,n);

% end of m-file