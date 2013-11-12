function hh=shade(tf,tr,trse,yd)
% SHADE  Plot shaded confidence bounds
%
% h=shade(t,fit,fitse,y)
%
% t: Time axis (t=[1:length(y)]').
% fit: Estimate (*).
% fitse: Standard error (*).
% y: Data (*).
%
% h: figure axis handle.
%
% Example: [fit,fitse]=dlr(y,z); shade([],fit,fitse,y)

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author for SHADE: Renata Romanowicz

% Graphing shell for plotting model estimates fit, shaded standard
% errors fitse and data y against time t, returning the graphics
% handle H.

if nargin==0
  disp(' ')
  disp(' SHADE  Plot shaded confidence bounds')
  disp(' ')
  disp(' h=shade(t,fit,fitse,y)')
  disp(' ')
  return
end

if nargin<1, tf=[]; end
if nargin<2, tr=[]; end
if nargin<3, trse=[]; end
if nargin<4, yd=[]; end

hh=shade0(tf,tr,trse,yd);

% end of m-file
