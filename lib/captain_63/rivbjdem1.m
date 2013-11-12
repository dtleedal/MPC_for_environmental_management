% RIVBJDEM1  Captain Toolbox demonstration
%
% Figures for supplementary documentation:
% "Upgraded TF Identification and Estimation Routines in CAPTAIN"
%
% Type 'rivbjdem1' at the command line to generate these figures.
%
% The script analyses the famous Box-Jenkins (BJ) gas furnace data
% using RIVBJ, RIVBJID, DTFM and DTFMOPT.

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact

set(0, 'defaulttextfontname', 'times', 'defaulttextfontsize', 12);
set(0, 'defaultaxesfontname', 'times', 'defaultaxesfontsize', 12);
set(0, 'defaultaxeslinewidth', 1)

disp(' ')
disp(' RIVBJDEM1  Captain Toolbox demonstration')
disp(' ')
disp(' This script analyses the famous Box-Jenkins (BJ) gas')
disp(' furnace data using RIVBJ, RIVBJID, DTFM and DTFMOPT.')
disp(' ')
disp(' It generates figures for the supplementary documentation:')
disp(' Upgraded TF Identification and Estimation Routines in CAPTAIN.')
disp(' ')
disp(' This demonstration script is intended for use with the')
disp(' above document and so is not annotated. However, you may')
disp(' press Ctrl + C to quit, followed by ''type rivbjdem1'' to examine')
disp(' the m-file, or ''edit rivbjdem1'' to open it for editing. In the')
disp(' latter case, it is possible to copy and paste lines into the')
disp(' Command Window as required.')
disp(' ')
disp(' --------------------------------------------------------')
disp('      Hit any key to continue or Ctrl + C to quit')
disp(' --------------------------------------------------------')
pause

load bjgas.dat
y=bjgas(:, 1)-mean(bjgas(:, 1));
u=bjgas(:, 2)-mean(bjgas(:, 2));

% automatically detect for SID Toolbox
% (or set sid=0 to test Captain without SID Toolbox)
sid=exist('iddemo.m')==2;

if sid  % using SID Toolbox

  disp(' ')
  disp('Using SID Toolbox for ARMA noise model estimation')
  disp(' ')

  D=iddata(y, u, 1);

  % Box-Jenkins (BJ) discrete-time (DT) model identification
  Md=rivbjid([y u], [1 1 1 2 0; 1 3 3 2 0], 1, [], [], sid)
  disp([Md.a(2) Md.b(4) Md.c(2:3); sqrt(diag(Md.CovarianceMatrix))'])
  ym=sim(D,Md);
  e=Md.UserData.e;
  eef=Md.UserData.er;
  figure(1)
  tab=acf(eef, 20);

  % [1 1 3 2 0] continuous-time (CT) model estimation
  %     (slightly better explanation of the data)
  Mc=rivcbj([y u], [1 1 3 2 0], [-1 1 1 0], -0.1, [], sid);
  disp([Mc.a(2) Mc.b Mc.c(2:3); sqrt(diag(Mc.CovarianceMatrix))'])
  ec=Mc.UserData.e;
  eefc=Mc.UserData.er;
  figure(2)
  tab=acf(eefc, 20);

  % time variable parameter (TVP) estimation of BJ model
  [nvr, opts]=dtfmopt(y, u, [1 3 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 3 3], 0, nvr);
  t=ones(size(y));
  t=cumsum(t);

  % Figure 3 from "Updated TF Identification and Estimation Routines"
  figure(3)
  plot(t,par)
  legend('a_1(k) estimate', 'b_0(k) estimate', 'b_1(k) estimate', 'b_2(k) estimate', 2)
  axis([min(t), max(t), -2, 2])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')

  % consider a shorter data set
  y=bjgas(1:160, 1)-mean(bjgas(1:160, 1));
  u=bjgas(1:160, 2)-mean(bjgas(1:160, 2));
  Ds=iddata(y, u, 1);

  % RIVBJID identification comparison
  Md=rivbjid([y u], [1 1 3 2 0; 1 3 3 2 1], 1, [], [], sid);
  disp([Md.a(2) Md.b(4) Md.c(2:3); sqrt(diag(Md.CovarianceMatrix))'])
  ym=sim(Ds, Md);
  e=y-ym.y;
  eef=filter(Md.c, Md.d, e);
  figure(4)
  tab=acf(eef, 20);

  % RIVBJ estimation of [1 1 3 2 0] model
  Md=rivbjid([y u], [1 1 3 2 0; 1 1 3 2 0], 1, [], [], sid);
  disp([Md.a(2) Md.b(4) Md.c(2:3); sqrt(diag(Md.CovarianceMatrix))'])
  ym=sim(Ds, Md);
  e=y-ym.y;
  eef=filter(Md.c, Md.d, e);
  figure(5)
  tab=acf(eef, 20);

  % TVP estimation of [1 1 3 2 0] model
  [nvr, opts]=dtfmopt(y, u, [1 1 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 1 3], 0, nvr);
  t=ones(size(y));
  t=cumsum(t);

  % Figure 4 from "Updated TF Identification and Estimation Routines"
  figure(6)
  subplot(121)
  plot(t, par)
  legend('a_1(k) estimate','b_0(k) estimate')
  axis([min(t), max(t), -1.5, 0])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')
  title('[1 1 3 0 0] Model')

  % TVP estimation of BJ model
  [nvr, opts]=dtfmopt(y, u, [1 3 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 3 3], 0, nvr);
  subplot(122)
  plot(t, par);
  legend('a_1(k) estimate', 'b_0(k) estimate', 'b_1(k) estimate', 'b_2(k) estimate', 4);
  axis([min(t), max(t), -1.5, 0])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')
  title('[1 3 3 0 0] Model')

  % RIVBJ estimation of [1 1 3 2 1] model (better residuals)
  Md=rivbjid([y u], [1 1 3 2 1; 1 1 3 2 1], 1, [], [], sid)
  disp('This is the final estimated constant parameter DT model in object code format')
  disp(' ')
  disp(' --------------------------------------------------------')
  disp('      Hit any key to continue or Ctrl + C to quit')
  disp(' --------------------------------------------------------')
  disp(' ')
  pause
  disp([Md.a(2) Md.b(4) Md.c(2:3) Md.d(2); sqrt(diag(Md.CovarianceMatrix))'])
  ym=sim(Ds, Md);
  e=y-ym.y;
  eef=filter(Md.c, Md.d, e);
  figure(7)
  tab=acf(eef,20);
  
else  % without using SID Toolbox

  % Box-Jenkins (BJ) discrete-time (DT) model identification
  [th, stats, e, RR, eef, var, Ps, Pn]=...
      rivbjid([y u], [1 1 1 2 0; 1 3 3 2 0], 1, [], [], sid);
  [ad, bd, cd, dd, Pd, td]=getparbj(th);
  disp([ad(2) bd(4) cd(2:3); sqrt(diag(Pd))'])
  figure(1)
  tab=acf(eef, 20);

  % [1 1 3 2 0] continuous-time (CT) model estimation
  [thc, statsc, ec, eefc, varc, Psc, Pnc]=...
      rivcbj([y u], [1 1 3 2 0], [-1 1 1 0], -0.1, [], sid);
  [ac, bc, cc, dc, Pc, tdc]=getparbj(thc);
  disp([ac(2) bc cc(2:3);sqrt(diag(Pc))'])
  figure(2)
  tab=acf(eefc, 20);

  % time variable parameter (TVP) estimation of BJ model
  [nvr, opts]=dtfmopt(y, u, [1 3 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 3 3], 0, nvr);
  t=ones(size(y));
  t=cumsum(t);

  % Figure 3 from "Updated TF Identification and Estimation Routines"
  figure(3)
  plot(t, par)
  legend('a_1(k) estimate','b_0(k) estimate','b_1(k) estimate','b_2(k) estimate', 2)
  axis([min(t), max(t), -2, 2])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')

  % consider a shorter data set
  y=bjgas(1:160, 1)-mean(bjgas(1:160, 1));
  u=bjgas(1:160, 2)-mean(bjgas(1:160, 2));

  % RIVBJID identification comparison
  [th, stats, e, RR, eef, var, Ps, Pn]=...
      rivbjid([y u], [1 1 3 2 0; 1 3 3 2 1], 1, [], [], sid);
  [ad, bd, cd, dd, Pd, td]=getparbj(th);
  disp([ad(2) bd(4) cd(2:3) dd(2); sqrt(diag(Pd))'])
  figure(4)
  tab=acf(eef, 20);

  % RIVBJ estimation of [1 1 3 2 0] model
  [th, stats, e, RR, eef, var, Ps, Pn]=rivbjid([y u], [1 1 3 2 0; 1 1 3 2 0], 1, [], [], sid);
  [ad, bd, cd, dd, Pd, td]=getparbj(th);
  disp([ad(2) bd(4) cd(2:3); sqrt(diag(Pd))'])
  figure(5)
  tab=acf(eef, 20);

  % TVP estimation of [1 1 3 2 0] model
  [nvr, opts]=dtfmopt(y, u, [1 1 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 1 3], 0, nvr);
  t=ones(size(y));
  t=cumsum(t);

  % Figure 4 from "Updated TF Identification and Estimation Routines"
  figure(6)
  subplot(121)
  plot(t, par);
  legend('a_1(k) estimate','b_0(k) estimate')
  axis([min(t), max(t), -1.5, 0])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')
  title('[1 1 3 0 0] Model')

  % TVP estimation of BJ model
  [nvr, opts]=dtfmopt(y, u, [1 3 3], 0, 'f1', -2);
  [tfs, fit, fitse, par, parse, e, y0]=dtfm(y, u, [1 3 3], 0, nvr);
  subplot(122)
  plot(t, par)
  legend('a_1(k) estimate','b_0(k) estimate','b_1(k) estimate','b_2(k) estimate', 4)
  axis([min(t), max(t), -1.5, 0])
  xlabel('Number of Samples')
  ylabel('TVP Estimate')
  title('[1 3 3 0 0] Model')

  % RIVBJ estimation of [1 1 3 2 1] model (better residuals)
  [th, stats, e, RR, eef, var, Ps, Pn]=...
      rivbjid([y u], [1 1 3 2 1; 1 1 3 2 1], 1, [], [], sid);
  [ad, bd, cd, dd, Pd, td]=getparbj(th);
  disp([ad(2) bd(4) cd(2:3) dd(2); sqrt(diag(Pd))'])
  figure(7)
  tab=acf(eef,20);

end

% end of m-file
