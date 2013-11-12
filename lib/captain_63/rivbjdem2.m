% RIVBJDEM2  Captain Toolbox demonstration
%
% Figures for supplementary documentation:
% "Upgraded TF Identification and Estimation Routines in CAPTAIN"
%
% Type 'rivbjdem2' at the command line to generate these figures.
%
% The script uses RIVBJ and RIVCBJ to show the advantages of
% continuous-time estimation in the case of rapidly sampled data
% and a stiff system

% Copyright (c) 2007 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

clear all
close all
format compact

set(0, 'defaulttextfontname', 'times', 'defaulttextfontsize', 12);
set(0, 'defaultaxesfontname', 'times', 'defaultaxesfontsize', 12);
set(0, 'defaultaxeslinewidth', 1)

disp(' ')
disp(' RIVBJDEM2  Captain Toolbox demonstration')
disp(' ')
disp(' This script uses RIVBJ and RIVCBJ to show the advantages of')
disp(' continuous-time (CT) estimation in the case of rapidly sampled')
disp(' data and a stiff system.')
disp(' ')
disp(' It generates figures for the supplementary documentation:')
disp(' Upgraded TF Identification and Estimation Routines in CAPTAIN.')
disp(' ')
disp(' This demonstration script is intended for use with the')
disp(' above document and so is not annotated. However, you may')
disp(' press Ctrl + C to quit, followed by ''type rivbjdem2'' to examine')
disp(' the m-file, or ''edit rivbjdem2'' to open it for editing. In the')
disp(' latter case, it is possible to copy and paste lines into the')
disp(' Command Window as required.')
disp(' ')
disp(' Warning: this demonstration is optimised for Matlab v7')
disp(' onwards and may yield poor results for earlier versions')
disp(' of Matlab.')
disp(' ')
disp(' --------------------------------------------------------')
disp('      Hit any key to continue or Ctrl + C to quit')
disp(' --------------------------------------------------------')
pause

load srivcbj.dat
y=srivcbj(1:5000,1);
u=srivcbj(1:5000,2);
x=srivcbj(1:5000,3);
n=srivcbj(1:5000,4);

Nc=[-120  -1560  3600];         % true numerator
Dc=[1  30.2  3607  750  3600];  % true denominator
Cd=[1 -1.4 0.7];                % true AR
Dd=[1 0.5];                     % true MA
Ts=0.005;                       % sampling interval

% automatically detect for SID Toolbox
% (or set sid=0 to test Captain without SID Toolbox)
sid=exist('iddemo.m')==2;

if sid  % using SID Toolbox

  disp(' ')
  disp('Using SID Toolbox for ARMA noise model estimation')
  disp(' ')

  U=iddata([], [u], 0.005);
  Mct=idpoly(Dc, Nc, [], [], [], [], 0) 

  % RIVCBJ estimation without noise model (SRIVC)
  Mc1=rivcbj([y u], [4 3 7 0 0], [-1 Ts 1 0], -0.1, [], sid);
  ac=Mc1.a;
  bc=Mc1.b;
  ym=sim(Mc1, U);
  Pc=Mc1.CovarianceMatrix;
  e=y-ym.y;
  disp(num2str([Dc(2:end) Nc; ac(2:end) bc; sqrt(diag(Pc))'])) 

  % IVARMA identification of ARMA noise model order from RIVCBJ generated e
  [p, q, C, D, Pcc, resvar, eef, RR]=ivarmaid(e, [1 1 3 3], 50, 1, sid);

  % Full RIVCBJ estimation using identified ARMA noise model
  Mc2=rivcbj([y u], [4 3 7 p q], [6 Ts 1 0], -0.1, [], sid)
  disp('This is the final estimated CT model in object code format')
  disp(' ')
  disp(' --------------------------------------------------------')
  disp('      Hit any key to continue or Ctrl + C to quit')
  disp(' --------------------------------------------------------')
  disp(' ')
  pause
  
  ac=Mc2.a;
  bc=Mc2.b;
  cc=Mc2.d;dc=Mc2.c;
  ym=y-Mc2.UserData.er;
  Pc=Mc2.CovarianceMatrix;
  e=y-ym;
  disp(num2str([Dc(2:end) Nc Cd(2:end) Dd(2:end); ac(2:end) bc cc(2:end) dc(2:end); sqrt(diag(Pc))']))
  % next two lines added to avoid noise model with incorrect CT
  % form causing instability in later step response calculation
  Mc2.c=1;
  Mc2.d=1;

  % Full discrete-time (DT) RIVBJ estimation and conversion to CT (poor results)
  Md=rivbj([y u], [4 4 7 p q], [0.005], [], [], [], sid);
  % Notes:
  % (a) This is not as good as the CT estimation because there is a pair of
  %     fairly slow roots and a pair of quite fast ones, so the slow ones
  %     will be oversampled (see roots(ad)); PEM does not estimate it at all,
  %     probably for this reason
  % (b) The orders of continuous and discrete TF don't have to be the
  %     same ([4 4 7 p q] is correct for DT)
  Md.Ts=0.005;  % since RIVBJ gives unity sampling interval Ts 
  ad=Md.a;
  bd=Md.b;
  cd=Md.d;
  dd=Md.c;
  ymd=sim(Md, U);
  Pc=Md.CovarianceMatrix;
  e=y-ymd.y;
  [bdc2, adc2]=d2cm(bd(8:11), ad, Ts, 'zoh');

  % PEM estimation and conversion to CT (poor results)
  % Note: [4 3 7] model used because other possible models are worse
  D=iddata(y, u, 0.005)
  TH=pem(D, [0 3 q p 4 7]);
  [bdc3, adc3]=d2cm(TH.b(8:10), TH.f, 0.005, 'zoh');

  % Figure 1 from "Updated TF Identification and Estimation Routines"
  figure(1)
  bode(Mct, 'r', Mc2, 'k', Md, 'b', TH, 'g')
  legend('True', 'RIVCBJ [4 3 7]', 'RIVBJ with D2CM [4 4 7]', 'PEM with D2CM [4 3 7]', 3)

  % Figure 2 from "Updated TF Identification and Estimation Routines"
  %   Complaints below from SID step function because of object coding
  %   differences in hybrid model (CT model with DT noise is not allowed in SID)
  figure(2)
  step(Mct, 'r', Mc2, 'k', Md, 'b');
  legend('True', 'RIVCBJ [4 3 7]', 'RIVBJ with D2CM [4 4 7]')

else  % without using SID Toolbox

  % RIVCBJ estimation without noise model (SRIVC)
  [th, stats, e, eef, var, Ps, Pn, y0]=...
      rivcbj([y u], [4 3 7 0 0], [-1 Ts 1 0], -0.1, [], sid);
  [ac, bc, cc, dc, Pc, delc]=getparbj(th);

  disp(num2str([Dc(2:end) Nc; ac(2:end) bc; sqrt(diag(Pc))'])) 

  %IVARMA identification of ARMA noise model order
  [p, q, C, D, Pcc, resvar, eef, RR]=ivarmaid(e, [1 1 3 3], 50, 1, sid);

  % Full RIVCBJ estimation
  [th, stats, e, eef, var, Ps, Pn, y0]=...
      rivcbj([y u], [4 3 7 p q], [6 Ts 1 0], -0.1, [], sid);
  [ac, bc, cc, dc, Pc, delc]=getparbj(th);
  disp(num2str([Dc(2:end) Nc Cd(2:end) Dd(2:end); ac(2:end) bc cc(2:end) dc(2:end); sqrt(diag(Pc))']))

  % Full RIVBJ estimation and conversion to CT (relatively poor results)
  [th, stats, e, eef, vr, Ps, Pn, y0]=rivbj([y u], [4 4 7 p q], [0.005], [], [], [], sid);
  % Note: see above object coded comment
  [ad, bd, cd, dc, Pd, deld]=getparbj(th);
  [bdc2, adc2] = md2cm(bd(8:11), ad, Ts);

  % Figure 1 from "Updated TF Identification and Estimation Routines"
  if exist('bode.m')==2
    figure(1)
    bode(Nc, Dc)
    hold on
    bode(bc, ac)
    bode(bdc2, adc2)
    legend('True', 'RIVCBJ [4 3 7]', 'RIVBJ with D2CM [4 4 7]');
  end

  % Figure 2 from "Updated TF Identification and Estimation Routines"
  if exist('step.m')==2
    figure(2)
    step(Nc, Dc, 'r')
    hold on
    step(bc, ac, 'k')
    step(bdc2, adc2, 'b')
    legend('True', 'RIVCBJ [4 3 7]', 'RIVBJ with D2CM [4 4 7]')
  end
end

% end of m-file
 

