function [A,B,C,D,L] = make_test_model(nL, dt,dz)

if (nargin == 0)
% define discretisation parameters
nL = 40;    % number of ocean layers
dt = 0.1;   % fraction of year for time step
dz = 100;   % depth of ocean layers
end
n_model=2; % return an MC ensemble of 2 members (then use first)
%Physical parameters
fnl_bl = 0.42;                      % fraction NH land baseline
fnl_dist = 1;                       % (1) uniform distribution
fnl_spread = 0;                    % +/- 10 percent range
fnl = IAGP_return_MC_par(n_model,fnl_bl,fnl_dist,fnl_spread);

fsl_bl = 0.21;                      % fraction SH land
fsl_dist = 1;                       % uniform
fsl_spread = 0;                    % +/- 10%
fsl = IAGP_return_MC_par(n_model,fsl_bl,fsl_dist,fsl_spread);


klo_bl = 1;                         % flux coefficient land/ocean (W/m^2/C)
klo_dist = 1;                       % uniform
klo_spread = 0;                    % +/- 10%
klo = IAGP_return_MC_par(n_model,klo_bl,klo_dist,klo_spread);


kns_bl = 1;                         % flux coefficient north/south (W/m^2/C)
kns_dist = 1;                       % uniform
kns_spread = 0;                    % +/- 10%
kns = IAGP_return_MC_par(n_model,kns_bl,kns_dist,kns_spread);


ll_bl = 0.716;                      % land longwave outgoing feedback parameter (W/m^2/C)
ll_dist = 1;                        % uniform
ll_spread = 30;                     % +/- 10%
ll = IAGP_return_MC_par(n_model,ll_bl,ll_dist,ll_spread);
%ll = [linspace(ll_bl*(100-ll_spread)/100,ll_bl*(100+2.*ll_spread)/100,5) linspace(ll_bl*(100-ll_spread)/100,ll_bl*(100+2.*ll_spread)/100,5)];

lo_bl = 3.648;                      % ocean longwave outgoing feedback parameter (W/m^2/C)
lo_dist = 1;                        % uniform
lo_spread = 30;                     % +/- 10%
lo = IAGP_return_MC_par(n_model,lo_bl,lo_dist,lo_spread);

%lo = [linspace(lo_bl*(100-lo_spread)/100,lo_bl*(100+2.*lo_spread)/100,5) linspace(lo_bl*(100-lo_spread)/100,lo_bl*(100+2.*lo_spread)/100,5)];

hm_bl = 90;                         % depth of mixed layer (m)
hm_dist = 1;                        % uniform
hm_spread = 0;                     % +/- 10%
hm = IAGP_return_MC_par(n_model,hm_bl,hm_dist,hm_spread);

rho = 1025.98;                      % density of sea water (kg/m3)
cp = 3989.8;                        % specific heat capacity of sea water (J/kg)

diffusivity_bl = 2.547e-2;           % diffusivity parameter
diffusivity_dist = 1;               % uniform
diffusivity_spread = 20;             % +/- 5%
diffusivity = IAGP_return_MC_par(n_model,diffusivity_bl,diffusivity_dist,diffusivity_spread);
%diffusivity = [linspace(diffusivity_bl*(100-diffusivity_spread)/100,diffusivity_bl*(100+diffusivity_spread)/100,5) linspace(diffusivity_bl*(100-diffusivity_spread)/100,diffusivity_bl*(100+diffusivity_spread)/100,6)];
%diffusivity = linspace(diffusivity_bl*(100-diffusivity_spread)/100,diffusivity_bl*(100+diffusivity_spread)/100,10);

kdo = (diffusivity).^2*3600*24*365*1;	% mixed layer to deep ocean diffusivity	(m2/y)

w_bl = 0;                           % overturning
w_dist = 1;                         % uniform
w_spread = 0;                       % +/- 0%
w = IAGP_return_MC_par(n_model,w_bl,w_dist,w_spread);

th_bl = 0;                          % upwelling
th_dist = 1;                        % uniform
th_spread = 0;                      % +/- 0%
th = IAGP_return_MC_par(n_model,th_bl,th_dist,th_spread);

% PI_bl = -2.9;                     % NH temp to sea ice extent coefficient
% PI_dist = 1;                        % uniform
% PI_spread = 0;                     % +/- 10%
% PI = IAGP_return_MC_par(n_model,PI_bl,PI_dist,PI_spread);
%
% phi_bl = -0.2;                      % SO2 to emissions to forcing coefficient
% phi_dist = 1;                       % uniform
% phi_spread = 0;                    % +/- 10%
% phi = IAGP_return_MC_par(n_model,phi_bl,phi_dist,phi_spread);

% NH ocean fraction (1 - northern hemisphere land)
fno = 1-fnl;

% SH ocean fraction (1 - southern hemisphere land)
fso = 1-fsl;

% effective bulk heat capacity of mixed layer (W*year/m^2/C)
cm = rho*cp.*hm./3600/24/365;

for mc = 2:n_model
    % build the state space model
    % build state space model
    %---------------------------------
    
    
    
    %Northern hemisphere
    AN = 1 + (klo(mc).^2/fno(mc)/cm(mc)/(klo(mc)+fnl(mc)*ll(mc))...
        - lo(mc)/cm(mc) - kdo(mc)/hm(mc)/0.5/dz...
        - w(mc)*th(mc)/hm(mc) - kns(mc)/fno(mc)/cm(mc) - klo(mc)/fno(mc)/cm(mc))*dt;
    BN = (1/cm(mc) + klo(mc)*fnl(mc)/fno(mc)/cm(mc)/(klo(mc)+fnl(mc)*ll(mc)))*dt;
    CN = (kdo(mc)/hm(mc)/0.5/dz + w(mc)/hm(mc))*dt;
    DN = (kns(mc)/cm(mc)/fno(mc))*dt;
    EN = klo(mc)/(fnl(mc)*ll(mc)+klo(mc));
    FN = fnl(mc)/(fnl(mc)*ll(mc)+klo(mc));
    
    %Southern hemisphere
    AS = 1 + (klo(mc)^2/fso(mc)/cm(mc)/(klo(mc)+fsl(mc)*ll(mc))...
        - lo(mc)/cm(mc) - kdo(mc)/hm(mc)/0.5/dz...
        - w(mc)*th(mc)/hm(mc) - kns(mc)/fso(mc)/cm(mc)...
        - klo(mc)/fso(mc)/cm(mc))*dt;
    BS = (1/cm(mc) + klo(mc)*fsl(mc)/fso(mc)/cm(mc)/(klo(mc)+fsl(mc)*ll(mc)))*dt;
    CS = (kdo(mc)/hm(mc)/0.5/dz + w(mc)/hm(mc))*dt;
    DS = (kns(mc)/cm(mc)/fso(mc))*dt;
    ES = klo(mc)/(fsl(mc)*ll(mc)+klo(mc));
    FS = fsl(mc)/(fsl(mc)*ll(mc)+klo(mc));
    
    %Ocean
    AO = (1-dt*(kdo(mc)/0.5/dz/dz+kdo(mc)/dz/dz+w(mc)/dz));
    BO = dt*(kdo(mc)/0.5/dz/dz);
    CO = dt*(kdo(mc)/dz/dz+w(mc)/dz);
    DO = (1-dt*(2*kdo(mc)/dz/dz+w(mc)/dz));
    EO = dt*(kdo(mc)/dz/dz);
    FO = dt*(kdo(mc)/dz/dz+w(mc)/dz);
    GO = (1-dt*(kdo(mc)/dz/dz+w(mc)/dz));
    HO = dt*(kdo(mc)/dz/dz);
    IO = dt*(th(mc)*w(mc)/dz);
    
    %------------------------------
    %State space version
    % construct the upper left block of the A matrix
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
    %Btop(1,3) = BN*1/phi(mc);
    % construct bottom half of B matrix
    Blower = zeros(nL,2);
    Blower(1,2) = BS;
    
    
    % put the B matrix together
    B = [Btop; Blower];
    
    
    % construct left side of C matrix together
    % the outputs are: NH layer 1, NH land, SH layer 1, SH land, NH mean, SH
    % mean, Global mean, sea ice extent
    
    % the original 8 state output C matrix
    %     Cleft = zeros(8,nL);
    %     Cleft(1,1) = 1;
    %     Cleft(2,1) = EN;
    %     Cleft(5,1) = (fno + fnl*EN);
    %     Cleft(7,1) = 0.5*(fno + fnl*EN);
    %     Cleft(8,1) = PI*fno;
    %     % construct right side of C matrix together
    %     Cright = zeros(8,nL);
    %     Cright(3,1) = 1;
    %     Cright(4,1) = ES;
    %     Cright(6,1) = (fso + fsl*ES);
    %     Cright(7,1) = 0.5*(fso + fsl*ES);
    %     % put together the C matrix
    %     C = [Cleft Cright];
    
    % the original 4 state C matrix
    %     Cleft = zeros(4,nL);
    %     Cleft(1,1) = (fno(mc) + fnl(mc)*EN);
    %     Cleft(3,1) = 0.5*(fno(mc) + fnl(mc)*EN);
    %     Cleft(4,1) = PI(mc)*(fno(mc) + fnl(mc)*EN);
    %     % construct right side of C matrix together
    %     Cright = zeros(4,nL);
    %     Cright(2,1) = (fso(mc) + fsl(mc)*ES);
    %     Cright(3,1) = 0.5*(fso(mc) + fsl(mc)*ES);
    %     % put together the C matrix
    %     C(:,:,mc) = [Cleft Cright];
    
    % the new 4 state C matrix for N S temp difference
    Cleft = zeros(2,nL);
    Cleft(1,1) = (fno(mc) + fnl(mc)*EN);
    
    % construct right side of C matrix together
    Cright = zeros(2,nL);
    Cright(2,1) = (fso(mc) + fsl(mc)*ES);
    
    % put together the C matrix
    C = [Cleft Cright];
    
    % build the D matrix
    % the old 8 state D matrix
    %     D = [0                   0                  0;...
    %         FN                   0                  FN*1/phi;...
    %         0                    0                  0;...
    %         0                    FS                 0;...
    %         (fnl*FN)         0                  (fnl*FN*1/phi);...
    %         0                    (fsl*FS)       0;...
    %         0.5*(fnl*FN)     0.5*(fsl*FS)   0.5*(fnl*FN*1/phi);...
    %         PI*(fnl*FN)  0                  PI*(fnl*FN*1/phi)];
    
    % the new 2 state D matrix
    D = [(fnl(mc)*FN)              0;...
        0   (fsl(mc)*FS)];
end
mc = 2;
start_year = 2;
start_year_less1 = start_year-1;
end_year = 5000;

t = (start_year_less1:end_year)';
forcing = zeros(size(t));
forcing(10:end) = 3.74;
all_X_states = zeros(2*nL,length(t));
all_Y_states = zeros(2,length(t));
for time_index = start_year:end_year
    
    
    last_state = all_X_states(:,time_index-1);
    
    
    % form forcing and emissions into a vector
    
    %    * * * * * * * * * * * * * * * * * * * * * * * * * *
    %    *                                                 *
    %    *                 **IMPORTANT**                   *
    %    *           THE FORCINGS ALL HAVE TO BE           *
    %    *                 MULTIPLIED BY 2                 *
    %    *           BECAUSE OF THE WAY MAGICC             *
    %    *            DEFINES THE HEMISPHERES              *
    %    *                                                 *
    %    * * * * * * * * * * * * * * * * * * * * * * * * * *
    forcings = [forcing(time_index); forcing(time_index)].*2;
    % do Kalman Filter step
    [X_k, Y_k] = stateSpaceSim(last_state,forcings,dt,A,B,C,D);
    
    
    all_X_states(:,time_index) = X_k;
    % place results into matrix
    all_Y_states(:,time_index) = Y_k;
    
    
    
end
L = mean(all_Y_states(:,end));
disp(['The LAMBDA is: ' num2str(mean(all_Y_states(:,end)))])
A = A(:,:,end);
B = B(:,:,end);
C = C(:,:,end);
D = D(:,:,end);
