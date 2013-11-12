function [y, meany, stdy] = stand0(x, meany, stdy)
% STAND  Standardise or de-standardise a matrix by columns
%
% [x,my,sy]=stand(y,my,sy)
%
% y: Time series (column matrix) (*)
% my: Vector of mean values for each column of y
% sy: Vector of standard deviations for each column of y
%
% x: Standardised (1 input argument) or de-standardised
%      (3 input arguments) equivalent of y
% my: Vector of mean values for each column of y
% sy: Vector of standard deviations for each column of y
%
% Examples:  A column vector y can be standardised or de-standardised
%   [x,my,sy]=stand(y) returns x=(y-my)./sy
%   [x]=stand(y,my,sy) returns x=y*sy+my
%
% See also FCAST, PREPZ

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% When the inputs my and sy are not supplied, this function 
% generates standardised variables by subtracting the mean from 
% each column of y and dividing by its standard deviation. In this 
% case, the function returns the standardised version of y in the 
% output argument x, together a vector of means (my) and 
% standard deviations (sy).
% 
% If, on the contrary, the input vectors/or scalars my and sy are 
% included in the call (i.e. there are three input arguments), then 
% the input matrix y is de-standardise column by column. In this 
% case, each column of y is multiplied by the standard deviation, 
% the mean is added and the de-standardised signal is returned as 
% the output argument x.

if nargin==0
  disp(' ')
  disp(' STAND  Standardise or de-standardise a matrix by columns')
  disp(' ')
  disp(' [x,my,sy]=stand(y,my,sy)')
  disp(' ')
  return
end

if nargin==1,
   [y,meany,stdy]=stand0(x);
elseif nargin==2, 
   [y,meany,stdy]=stand0(x,meany);
elseif nargin==3
   [y,meany,stdy]=stand0(x,meany,stdy);
end

% end of m-file