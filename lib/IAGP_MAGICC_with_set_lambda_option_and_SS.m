%MAGICCs test
%A.Jarvis, D.Leedal - 2013-14

% Input forcing
% 3.7Wm^-2 is the generally accepted forcing resulting from
% a doubling of atmospheric CO2 from pre-industrial levels
%forcing = zeros(5000,1);forcing(21:end) = 3.7;
forcing = zeros(100,1);
pre_industrial_CO2 = 278;
CO2_conc = pre_industrial_CO2;

% make the forcing series for a 1% per year CO2 increase to pre-ind doubling 
for i = 1:100
    forcing(i) = 5.35*log(CO2_conc/278);
    CO2_conc = CO2_conc + 0.01*CO2_conc;
    if (CO2_conc > 2*pre_industrial_CO2); CO2_conc = 2*pre_industrial_CO2; end
end
doubling_year = find(forcing == 5.35*log((2*pre_industrial_CO2)/pre_industrial_CO2));
forcing = forcing(1:doubling_year(1));
%---------------------------------
%Physical parameters
fnl = (0.42.*1)/1;                      % fraction Norther Hemisphere (NH) land
fno = (1-fnl)/1;                        % fraction of NH ocean
fsl = (0.21.*1)/1;                      % fraction of SH land
fso = (1-fsl)/1;                        % fraction of SH ocean
klo = 1;                                % flux coefficient land/ocean (W/m^2/C)
kns = 1;                                % flux coefficient north/south (W/m^2/C)
% ll and lo are now set in the routine at the top
% but the previous values are kept here for comparison
%ll = 0.716;                            % land longwave outgoing feedback parameter (W/m^2/C)
%lo = 3.648;                            % ocean longwave outgoing feedback parameter (W/m^2/C)
hm = 90;                                % depth of mixed layer (m)
rho = 1025.98; 							%density of sea water (kg/m^3)
cp = 3989.8;							%specific heat capacity of sea water (J/kg)
cm = (rho*cp*hm/3600/24/365);           %specific heat capacity of the ocean
kdo = (2.547*1e-2)^2*3600*24*365*1;     % mixed layer to deep ocean diffusivity	(m^2/y)
%kdo = 3155.76;                         % this is the value used in the RIVM report
kdo = 3155.76*1.5;
w = 4;                                  % overturning (not used)
th = 1;                                 % upwelling (not used)

%---------------------------------
% Select the required 2xCO2 equilibrium climate sensitivity and
% the ratio of equilibrium land to ocean temperature (RLO).
% Not all combinations of values can be achieved:
% if ll (lambda land) comes back negative
% then choose an alternative 2xCO2 equilibrium climate sensitivity and/or RLO
eq_climate_sensitivity = 2.5; % 2xCO2 Eq. climate sensitivity (K)
lambda = 3.7/eq_climate_sensitivity;% In TF terms, lambda is the SSG
RLO = 1.2; % = (eq. temp. land)/(eq. temp. ocean)
fo = (fno + fso)/1; % fraction of ocean
fl = (fnl + fsl)/1; % fracion of land
%global PARS % set up a global structure for the optimization function to access
PARS.lambda = lambda;
PARS.fno = fno;
PARS.fso = fso;
PARS.fnl = fnl;
PARS.fsl = fsl;
PARS.fo = fo;
PARS.fl = fl;
PARS.RLO = RLO;
PARS.klo = klo;
PARS.kns = kns;
% call fminsearchbnd with optimize_land_and_ocean_lambdas to set the value
% of lo that generates the best value of ll in order to
% make the RLO = the chosen value
lo = fminsearch('IAGP_optimize_land_and_ocean_lambdas',0.8,[],PARS);
% repeat the body of optimize_land_and_ocean_lambdas
% using the returned value of lo as input to get the final
% ll and to check the temperature ratio
ll = (lambda*(fo + fl*RLO) - fo*lo)/(fl*RLO);
C_sens = [...
    (fno*lo + klo + kns)    -klo            -kns                0;
    -klo                    (fnl*ll + klo)  0                   0;
    -kns                    0               fso*lo + klo + kns  -klo
    0                       0               -klo                (fsl*ll + klo)];
Q_sens = [fno fnl fso fsl]'*3.7;
T_sens = C_sens\Q_sens;
t_land = (fnl*T_sens(2) + fsl*T_sens(4))/(fnl+fsl);
t_ocean = (fno*T_sens(1) + fso*T_sens(3))/(fno+fso);
disp(['lo = ' num2str(lo)])
disp(['ll = ' num2str(ll)])
disp(['RLO = ' num2str(t_land/t_ocean)]);

PARS2.forcing = forcing;
PARS2.lo = lo;
PARS2.ll = ll;
PARS2.TCR = 1.25;
kdo = fminsearch('IAGP_optimize_downwelling_for_TCR',3155,[],PARS2)
%---------------------------------
%Discretisation parameters
nL = 40;    % number of ocean layers
dt = 0.1;   % year fraction for time stepping
dz = 100;   % depth of each ocean layer in meters

%---------------------------------
%Aggregate parameters

%Northern hemisphere
AN = 1 + (klo^2/fno/cm/(klo+fnl*ll) - lo/cm - kdo/hm/0.5/dz - w*th/hm - kns/fno/cm - klo/fno/cm)*dt;
BN = (1/cm + klo*fnl/fno/cm/(klo+fnl*ll))*dt;
CN = (kdo/hm/0.5/dz + w/hm)*dt;
DN = (kns/cm/fno)*dt;
EN = klo/(fnl*ll+klo);
FN = fnl/(fnl*ll+klo);

%Southern hemisphere
AS = 1 + (klo^2/fso/cm/(klo+fsl*ll) - lo/cm - kdo/hm/0.5/dz - w*th/hm - kns/fso/cm - klo/fso/cm)*dt;
BS = (1/cm + klo*fsl/fso/cm/(klo+fsl*ll))*dt;
CS = (kdo/hm/0.5/dz + w/hm)*dt;
DS = (kns/cm/fso)*dt;
ES = klo/(fsl*ll+klo);
FS = fsl/(fsl*ll+klo);

%Ocean
%AO = (1-dt*(kdo/0.5/dz/dz+kdo/dz/dz+w/dz));
AO = (1-dt*(kdo/(0.5*dz*dz)+kdo/(dz*dz)+w/dz));
BO = dt*(kdo/0.5/dz/dz);
CO = dt*(kdo/dz/dz+w/dz);
DO = (1-dt*(2*kdo/dz/dz+w/dz));
EO = dt*(kdo/dz/dz);
FO = dt*(kdo/dz/dz+w/dz);
GO = (1-dt*(kdo/dz/dz+w/dz));
HO = dt*(kdo/dz/dz);
IO = dt*(th*w/dz);

%---------------------------------
%Re-sample disturbance
% QN = interp(QN,1/dt);
% QS = interp(QS,1/dt);
forcing = zeros(5000,1);forcing(21:end) = 3.7;
QN = reshape((forcing*ones(1,1/dt))',length(forcing).*(1/dt),1);
QS = QN;
t = (1:length(QN))'*dt; % make a time base in years
%---------------------------------
%Variables
TNO = zeros(length(t),nL); % temperature of the Northern Hemisphere ocean layers
TNL = zeros(size(t));      % temperature of the Northern Hemisphere land
TSO = zeros(length(t),nL); % temperature of the Southern Hemisphere ocean layers
TSL = zeros(size(t));      % temperature of the Southern Hemisphere land

%---------------------------------
%Simulate

for i = 2:length(t)
    
    for j = 1:nL
        
        if j==1
            
            %Northern hemisphere mixed layer
            TNO(i,j) = AN*TNO(i-1) + BN*QN(i) + CN*TNO(i-1,j+1) + DN*TSO(i-1,j);
            
            %Northern hemisphere land
            TNL(i) = EN*TNO(i) + FN*QN(i);
            
            %Southern hemisphere mixed layer
            TSO(i,j) = AS*TSO(i-1) + BS*QS(i) + CS*TSO(i-1,j+1) + DS*TNO(i-1,j);
            
            %Southern hemisphere land
            TSL(i) = ES*TSO(i) + FS*QS(i);
            
        elseif j==2
            
            TNO(i,j) = AO*TNO(i-1,j) + BO*TNO(i-1,j-1) + CO*TNO(i-1,j+1);
            
            TSO(i,j) = AO*TSO(i-1,j) + BO*TSO(i-1,j-1) + CO*TSO(i-1,j+1);
            
        elseif j>2 && j<nL
            
            TNO(i,j) = DO*TNO(i-1,j) + EO*TNO(i-1,j-1) + FO*TNO(i-1,j+1);
            
            TSO(i,j) = DO*TSO(i-1,j) + EO*TSO(i-1,j-1) + FO*TSO(i-1,j+1);
            
        elseif j==nL
            
            TNO(i,j) = GO*TNO(i-1,j) + HO*TNO(i-1,j-1) + IO*TNO(i-1,1);
            
            TSO(i,j) = GO*TSO(i-1,j) + HO*TSO(i-1,j-1) + IO*TSO(i-1,1);
            
        end
        
    end
    
end

%------------------------------
%Hemispheric averaging
TNh = fnl*TNL + fno*TNO(:,1);   % averaging by fraction of land and ocean
TSh = fsl*TSL + fso*TSO(:,1);   % same for SH

%------------------------------
%Subsampling
TNh = TNh(1:1/dt:end);          % put back to annual rather than 0.1 year samples
TSh = TSh(1:1/dt:end);
TGh = (TNh + TSh)/2;            % global temperature average
tLabel = t(1:1/dt:end);

%------------------------------
% plot the scalar model output
figure(3)
plot(tLabel,[TNh TSh TGh])
legend('NH surface temp. mean','SH surface temp. mean','global temp. mean')
title('scalar model structure')

%------------------------------
% move on to the equivalent state space representation
% construct the new A matrix incorporating 2 extra states to allow for
% effects of temp diff between N and S hemisphere
Aul = zeros(nL,nL);
Aul(1,[1 2]) = [AN CN];
Aul(end,1) = IO;
Aul(end,[end-1 end]) = [HO GO];
Aul(2,[1 2 3]) = [BO AO CO];
for i = 3:nL-1
    Aul(i,[i-1 i i+1]) = [EO DO FO];
end
% construct the lower right block of the A matrix
Alr = zeros(nL,nL);
Alr(1,[1 2]) = [AS CS];
Alr(end,1) = IO;
Alr(end,[end-1 end]) = [HO GO];
Alr(2,[1 2 3]) = [BO AO CO];
for i = 3:nL-1
    Alr(i,[i-1 i i+1]) = [EO DO FO];
end
% Construct the upper right blocks of the A matrix
Aur = zeros(nL,nL);
Aur(1,1) = DN;
% Construct the lower left blocks of the A matrix
All = zeros(nL,nL);
All(1,1) = DS;

% put it all together
A = [Aul Aur;All Alr];

% the new state space MagicC does not have the SO2 as a separate input
% or the sea ice as a separate output
% construct top half of B matrix
Btop = zeros(nL,2);
Btop(1,1) = BN;

% construct bottom half of B matrix
Blower = zeros(nL,2);
Blower(1,2) = BS;


% put the B matrix together
B = [Btop; Blower];



Cleft = zeros(2,nL);
Cleft(1,1) = (fno + fnl*EN);

% here the NH, SH averaging is done by the C and D matrices

% construct right side of C matrix together
Cright = zeros(2,nL);
Cright(2,1) = (fso + fsl*ES);

% put together the C matrix
C = [Cleft Cright];


% the new 2 state D matrix
D = [(fnl*FN)              0;...
    0   (fsl*FS)];

% allocate space for results
X = zeros(2*nL,length(t));  % this will hold the 40 ocean layer for each hemisphere 
Y = zeros(2,length(t));     % this will hold the mean NH and SH temperatures

% simulate
for i = 2:length(t)
    X(:,i) = A*X(:,i-1) + B*[QN(i);QS(i)];
    Y(:,i) = C*X(:,i) + D*[QN(i);QS(i)];
end


%------------------------------
%Subsampling for the state space output
TNhSS = Y(1,1:1/dt:end)';
TShSS = Y(2,1:1/dt:end)';
TGhSS = (TNhSS + TShSS)/2;


%------------------------------
% plot the state space output (should be the same as the scalar version!)
figure(3)
plot(tLabel,[TNhSS TShSS TGhSS])
legend('NH surface temp. mean','SH surface temp. mean','global temp. mean')
title('state space model structure')

%------------------------------
% look at the ocean layers for fun
figure(4)
plot(X(1:40,:)','k')
hold on
plot(X(41:80,:)','r')
title('NH (black) and SH (red) ocean layers') 




