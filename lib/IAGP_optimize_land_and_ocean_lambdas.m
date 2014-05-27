function err = IAGP_optimize_land_and_ocean_lambdas(lo,PARS)
%global PARS
lambda = PARS.lambda;
fno = PARS.fno;
fso = PARS.fso;
fnl = PARS.fnl;
fsl = PARS.fsl;
RLO = PARS.RLO;
fo = PARS.fo;
fl = PARS.fl;
klo = PARS.klo;
kns = PARS.kns;

ll = (lambda*(fo + fl*RLO) - fo*lo)/(fl*RLO);
C_sens = [...
    (fno*lo + klo + kns)    -klo            -kns                0;
    -klo                    (fnl*ll + klo)  0                   0;
    -kns                    0               (fso*lo + klo + kns)  -klo
    0                       0               -klo                (fsl*ll + klo)];
Q_sens = [fno fnl fso fsl]'*3.7;
T_sens = C_sens\Q_sens;
t_land = (fnl*T_sens(2) + fsl*T_sens(4))/(fnl+fsl);
t_ocean = (fno*T_sens(1) + fso*T_sens(3))/(fno+fso);
err = (t_land/t_ocean - RLO).^2;

end

