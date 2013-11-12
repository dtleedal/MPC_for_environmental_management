function [acl,bcl,bclu,bclv]=pipcl(a,b,v,am,bm)
% PIPCL  Closed loop transfer functions for PIP control
%
% [acl,bcl,bclu,bclv]=pipcl(a,b,v,am,bm)
%
% b/a: System polynomials (trunciated form) (*)
% v: PIP controller coefficients (see 'help gains') (*)
% bm/am: if bm/am are omitted, the PIP feedback structure
%          is used, otherwise these are the control model
%          polynomials (truncated form) required for the
%          PIP forward path structure with mismatch
%
% bcl/acl: Command to output (Yd => y)
% bclu/acl: Command to input (Yd => u)
% bclv/acl: Load disturbance, y=(bcl/acl).Yd+(bclv/acl).V
%
% See also PIP, PIPOPT, PIPCOM, PIPLIB, GAINS

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% For details of the feedback and forward path control structures
% see e.g. Taylor, C.J., Chotai, A. and Young, P.C., (1998),
% Proportional-Integral-Plus (PIP) control of time delay systems,
% Proceedings of the Institution of Mechanical Engineers,
% Journal of Systems and Control Engineering, 212, Part I, 37-48.
%
% In the brief help message above, Yd represents the command
% input, y the controlled output and u the control input, while
% bcl, acl, bclu and bclv are various polynomials in the backward
% shift operator. These can be used with the Matlab FILTER function
% to find the closed-loop input and output variables, i.e. without
% disturbances the closed-loop transfer functions are as follows:
%
% u = (bclu/acl).Yd
% y = (bcl/acl).Yd
%
% For a load disturbance V, the output is:
%
% y = (bcl/acl).Yd + (bclv/acl).V
%
% To investigate other types of disturbances, multivariable
% systems or the affect of disturbances on the input signal,
% use the Captain Toolbox Simulink library: PIPLIB.
%
% Note that {am,bm} represent the control model and so are not
% required for the feedback PIP control structure. In the case
% of the forward path structure, the internal model is given
% by the transfer function bm/am.
%
% Finally, the output arguments {acl,bcl,bclu,bclv} are all
% polynomials represented in their "full" form, including the
% leading unity for the denominator and unit time delay for
% the numerator. By contrast, the input polynomials {a,b,am,bm}
% should all be supplied in "truncated form", i.e. not including
% the leading unity and with one less time delay. This follows
% the notation of other PIP functions such as PIPOPT.

if nargin==0
  disp(' ')
  disp(' PIPCL  Closed loop transfer functions for PIP control')
  disp(' ')
  disp(' [acl,bcl,bclu,bclv]=pipcl(a,b,v,am,bm)')
  disp(' ')
  return
end

if nargin<1, a=[]; end
if nargin<2, b=[]; end
if nargin<3, v=[]; end
if nargin<4, am=[]; end
if nargin<5, bm=[]; end

if nargin<4
  [acl,bcl,bclu,bclv]=pipcl0(a,b,v);
else
  [acl,bcl,bclu,bclv]=pipcl0(a,b,v,am,bm);
end

% end of m-file