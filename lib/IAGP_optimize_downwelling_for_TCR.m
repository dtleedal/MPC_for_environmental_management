function [ err ] = IAGP_optimize_downwelling_for_TCR( par_set, PARS )
kdo = par_set(1);
w = par_set(2);
fnl = (0.42.*1)/1;                      % fraction Norther Hemisphere (NH) land
fno = (1-fnl)/1;                        % fraction of NH ocean
fsl = (0.21.*1)/1;                      % fraction of SH land
fso = (1-fsl)/1;                        % fraction of SH ocean
klo = 1;                                % flux coefficient land/ocean (W/m^2/C)
kns = 1;                                % flux coefficient north/south (W/m^2/C)
% ll and lo are now set in the routine at the top
% but the previous values are kept here for comparison
ll = PARS.ll;                            % land longwave outgoing feedback parameter (W/m^2/C)
lo = PARS.lo;                            % ocean longwave outgoing feedback parameter (W/m^2/C)
hm = 90;                                % depth of mixed layer (m)
rho = 1025.98; 							%density of sea water (kg/m^3)
cp = 3989.8;							%specific heat capacity of sea water (J/kg)
cm = (rho*cp*hm/3600/24/365);           %specific heat capacity of the ocean
%kdo is passed into the function        % mixed layer to deep ocean diffusivity	(m^2/y)

%w = 8;                                  % overturning (not used)
th = 1;                                 % upwelling (not used)



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
forcing = PARS.forcing;
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
err = (TGh(end) - PARS.TCR).^2;


end

