function Param = IAGP_model_setup(base_directory,project_name)
% Param = IAGP_model_setup(base_directory,project_name)
% base_directory    path to the home MPC directory that contains the lib,
%                   initializations, and test subdirectories
% project_name      the name of the project that will use this parameter set
%                   This is a directory within base_directory and was setup
%                   by a call to the IAGP_directory_setup function
% Param             A structure containing all the variables
%script to set up the controller parameters
% 
% D.T. Leedal
% 3-Sep-2013
%
% 
%

disp(['project name: '  project_name])

%% Load the spin up data to get the start date and check spin up exists
try
    load([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'spin_up_data.mat'])
    spin_up_start_year = spin_up_data(1,1);
    last_spinup_year = spin_up_data(end,1);
catch
    err_msg = 'File load failed. Can''t load the spin up data file. There should\nbe a .mat file called ''spin_up_data.mat'' containing a matrix called ''spin_up_data''\nsaved in the project''s ''inputs_and_outputs'' directory. Please check.';
    error('errorMsg:converted',err_msg)
end


%% sets project name variable
Param.project_name = project_name;

%% set make plots flag
Param.use_plotting_flag = 0;

%% set use the diary function flag
% set this to 1 if you want keep a text file of all the screen output
Param.use_diary_flag = 0;

%% use performance weighting flag
Param.use_weighting = 1;
%% define initial estimate of NH temp to sea ice and emissions to forcing coefficients
Param.PI_bl = -2.1;  % (-2.1)  temp --> sea ice
Param.phi_bl = -0.16; % (-0.16) emissions to forcing

%% set fitted DI_des parameters and flag
% if the use_fitted_DI_des_flag is set to 0 then there must be a target
% sea ice trajectory data file called sea_ice_target saved in file sea_ice_target
% the file must contain a column vector of sea ice minimum extents from the first year the control
% begins to the last year of the experiment + forecast horizon. The extents
% should be absolute not perturbations as the sea_ice_baseline will be
% subtracted from the values
%
% If the flag is set to 1, the IAGPperformMPC fits a smooth spline from the
% preceeding sea ice trend and stabilizing at the (absolute) target value by the target year
Param.use_fitted_DI_des_flag = 1;
Param.ice_stabilisation_level = 5;
Param.ice_stabilisation_year = 2040;

%% define GCM NH and SH temp. baseline and sea ice baseline
% set these as near as possible to the unperturbed (pre GHG) values for your model 
Param.NH_baseline = 13.111525367110948; % (13.1115 for HadGEM2)
Param.SH_baseline = 13.547842406021266; % (13.5478 for HadGEM2)
Param.sea_ice_baseline = 5.5;%pre-industrial minimum annaul sea ice extent (millions sqr. km) (5.5 for HadGEM2)
Param.global_baseline = mean([Param.NH_baseline Param.SH_baseline]); % (13.3297 for HadGEM2)

%% integral of error wind up limit (absolute value)
Param.integral_wind_up_limit = 50; % (50)

%% define the optimization weights
Param.optimize_weights.w1 = 100; % error between ice and target (100)
Param.optimize_weights.w2 = 1; % smoothness of ice (1)
Param.optimize_weights.w3 = 100; % integral of error (100)
Param.optimize_weights.w4 = 0; % target overshoot changed from zero in 2044 (0)
Param.optimize_weights.w5 = 0.05; % smoothness of emissions (0.05)
%% Set whether to use the present observation in the optimization
% if this is set to 0 then it's possible to do a true model inversion run
Param.use_last_observation = 1;
%% Define size of model ensemble and the emergent properties of the models
% this must be an even number and be at least 4 preferably around 10

% This section sets the 2 x CO2 eq. clim. sens. and trans. clim. resp. for the control model
% ensemble. The values are used later when specifying the feedback parameters for land and ocean, the ocean
% diffusivity, the downwelling, and the ratio of land to ocean equilibrium temperature (RLO)  parameters for MagicC.
% The functions attempt to find values for these 4 MagicC parameters that will achieve the specfied
% 2 x CO2 eq. clim. sens. and trans. clim. resp. for each of the control models
Param.n_model = 10;
Param.two_x_CO2_eq_T = repmat(linspace(2.5,3.5,Param.n_model/2),1,2);
Param.trans_clim_response = [Param.two_x_CO2_eq_T(1:Param.n_model/2).*0.5 Param.two_x_CO2_eq_T(Param.n_model/2+1:end).*0.6];
Param.RLO = [1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2];
if (length(Param.two_x_CO2_eq_T) ~= Param.n_model) || (length(Param.trans_clim_response) ~= Param.n_model)
    error('number of models does not match number of eq. sens. or tcr values')
end

%% define length of forecast horizon
% MPC uses a prediction horizon at each year step looking forward N steps
% somewhere between 5 and 10 works best here
Param.N = 6; % (5)

%% the year when controller kicks in 
% this is one more than the last year of the spin_up_data file 
Param.control_start = last_spinup_year+1;

%% the sum of squared errors weighting to emphasise recent fit
% i.e., S(t) = SSE_weight * S(t-1) + e(t)^2
Param.SSE_weight = 0.9;

%% the exponential weighting to apply to the sea ice fit
% to emphasise recent fit
Param.RMSD_exponential_weight = 0.92;

%% define simulation years
% start year is read from the spin_up_data file
Param.start_year = spin_up_start_year;
start_year_less1 = Param.start_year-1;
Param.end_year = 2100;
end_year_plusN = Param.end_year+Param.N;
Param.t = (start_year_less1:end_year_plusN)';

%% forcing conversion function parameters
% these are of the form:
% AerosolAndCloud_forcing = forcing_conversion(1)*GHG_forcing +
% forcing_conversion(2)
% simultanious equations from Andy:
% -0.2 = m*0.2 + c (1) gives zero forcing in 1860 for GHG offset of 0.2
% -1.3 = m*2.9 + c (2) derived from Faero (Forster et al 2007)
% m = -0.41
% c = -0.12

% you will now need to switch this in year 2037 to move from the original
% conversion function to the new one.
%% set the year when TOA forcing etc become available to calculate Faero
Param.Faero_instruments_begin = 1995;

%% set a switch for what to do with non-ghg anthropogenic forcing after instrument record begins
% (1) continue to use linear relation to ghg forcing
% (2) use estimate based on instrumental data and observed temperature
% (3) use a pre-specified value
Param.which_f_aero_to_use = 2;
% set value to use if option 3 above
Param.f_aero_pre_specified = -1.1; % -0.5 is a sensible value

%% prior to Faero_instruments_begin we use a functional estimate parameterised with
Param.forcing_conversion = [-0.41 -0.12];% 

%% set up initial outgoing SW radiation value
Param.SW_out_initial = 99.035; % (99.035) W/m^2 used to calculate albedo-induced forcing perturbation

%% the NH SH difference function
% NH_SH_difference = [-1.7282 -0.9875 -0.2071]; in development

%% now define schedule and space for node values for time-varying PI_bl
% first year should be first year of simulation
Param.use_PI_updating = 0; % keep switched off - in development!

%% set flag for PID-only control
Param.PID_only = 0; % set to 1 to do a simple PID run
Param.P_weight = 1;
Param.I_weight = 1;
Param.D_weight = 1;
%% adaptive gain parameters
% switch for temperature adaptive gain
Param.use_T_adaptive_gain_flag = 1;

% switch for sea ice adaptive gain
Param.use_SI_adaptive_gain_flag = 1;

% NVR for temperature adaptive gain
Param.a_g_NVR_T = 0.00000001; % (0.00005) the larger the more quickly the parameter can change

% NVR for sea ice adaptive gain
Param.a_g_NVR_SI = 0.0001; % (0.00001)

% initial adaptive gain
Param.a_g_0 = 1; % (1)

% initial error variance estimate
Param.P_k_0 = 10; % (10)

% set the adaptive gain dummy offset
% adaptive gain doesn't work close to zero so always add on this then subtract 
Param.offset = 10; % (10) 

%% Specify how volcanoes are treated
% The volcanic input can be in one of two forms:
% (1) a Tg value of stratospheric injection as used by Leeds
% (2) a volcanic dust index as used by Ben Kravitz
% Both methods have an equal residence time (e-folding time) in years

% switch for which type of volcano input
% 1 = Tg stratospheric injection
% 2 = volcanic dust index
% ~=(1 or 2) sets all volcanic forcing to zero
Param.use_volcano = 1; 
Param.volcanic_residence_time = 1; % (1)

if Param.use_volcano == 1
    % volcano input specifies a Tg of stratospheric input
    % so we use the same forcing coef as for geoengineering SO2
    Param.volcanic_forcing_coef = Param.phi_bl;
elseif Param.use_volcano == 2
    % volcano input specifies a dust index
    Param.volcanic_forcing_coef = -24;
else
    % volcano input is not used
    Param.volcanic_forcing_coef = 0;
end



%% define plant failure via weibul parameters
Param.use_plant_failure = 0; % flag to use or by-pass plant failure mode 
Param.weibul_a = 20;    % (20)
Param.weibul_b = 2;     % (2)

%% choose whether to make file backups (faster without)
Param.do_backup = 0;

%% Fixed physical parameters for MagicC
fnl = 0.42;                         % fraction NH land baseline
fno = 1-fnl;                        % NH ocean fraction (1 - northern hemisphere land)
fsl = 0.21;                         % fraction SH land
fso = 1-fsl;                        % SH ocean fraction (1 - southern hemisphere land)
klo = 1;                            % flux coefficient land/ocean (W/m^2/C)
kns = 1;                            % flux coefficient north/south (W/m^2/C)

hm = 90;                            % depth of mixed layer (m)

rho = 1025.98;                      % density of sea water (kg/m3)
cp = 3989.8;                        % specific heat capacity of sea water (J/kg)
cm = rho*cp.*hm./3600/24/365;       % effective bulk heat capacity of mixed layer (W*year/m^2/C)

th = 1;                             % upwelling


%% define MagicC discretisation parameters
Param.nL = 40;    % number of ocean layers
Param.dt = 0.1;   % fraction of year for time step
Param.dz = 100;   % depth of ocean layers 
%% Build state space models

dz = Param.dz;
nL = Param.nL;
dt = Param.dt;

for mc = 1:Param.n_model
    %% Ensemble physical parameters for MagicC
    
    VARIABLE_PARS = IAGP_MAGICC_get_lambdas_and_diffusivity(Param.two_x_CO2_eq_T(mc),Param.trans_clim_response(mc),0,Param.RLO(mc)) 
    ll = VARIABLE_PARS.ll;              % land longwave outgoing feedback parameter (W/m^2/C)
    lo = VARIABLE_PARS.lo;              % ocean longwave outgoing feedback parameter (W/m^2/C)
    w = VARIABLE_PARS.w;                % overturning
    kdo = VARIABLE_PARS.kdo;            % mixed layer to deep ocean diffusivity	(m^2/y)
    % build the state space model
    % build state space model
    %---------------------------------
    
    
    %Northern hemisphere
    AN = 1 + (klo.^2/fno/cm/(klo+fnl*ll)...
        - lo/cm - kdo/hm/0.5/dz...
        - w*th/hm - kns/fno/cm - klo/fno/cm)*dt;
    BN = (1/cm + klo*fnl/fno/cm/(klo+fnl*ll))*dt;
    CN = (kdo/hm/0.5/dz + w/hm)*dt;
    DN = (kns/cm/fno)*dt;
    EN = klo/(fnl*ll+klo);
    FN = fnl/(fnl*ll+klo);
    
    %Southern hemisphere
    AS = 1 + (klo^2/fso/cm/(klo+fsl*ll)...
        - lo/cm - kdo/hm/0.5/dz...
        - w*th/hm - kns/fso/cm...
        - klo/fso/cm)*dt;
    BS = (1/cm + klo*fsl/fso/cm/(klo+fsl*ll))*dt;
    CS = (kdo/hm/0.5/dz + w/hm)*dt;
    DS = (kns/cm/fso)*dt;
    ES = klo/(fsl*ll+klo);
    FS = fsl/(fsl*ll+klo);
    
    %Ocean
    AO = (1-dt*(kdo/0.5/dz/dz+kdo/dz/dz+w/dz));
    BO = dt*(kdo/0.5/dz/dz);
    CO = dt*(kdo/dz/dz+w/dz);
    DO = (1-dt*(2*kdo/dz/dz+w/dz));
    EO = dt*(kdo/dz/dz);
    FO = dt*(kdo/dz/dz+w/dz);
    GO = (1-dt*(kdo/dz/dz+w/dz));
    HO = dt*(kdo/dz/dz);
    IO = dt*(th*w/dz);
    
    %------------------------------
    %State space version
    % construct the upper left block of the A matrix
    
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
    Param.A(:,:,mc) = [Aul Aur;All Alr];
    
    % the new state space MagicC does not have the SO2 as a separate input
    % or the sea ice as a separate output
    % construct top half of B matrix
    Btop = zeros(nL,2);
    Btop(1,1) = BN;
    
    % construct bottom half of B matrix
    Blower = zeros(nL,2);
    Blower(1,2) = BS;
    
    
    % put the B matrix together
    Param.B(:,:,mc) = [Btop; Blower];
    
    
    
    Cleft = zeros(2,nL);
    Cleft(1,1) = (fno + fnl*EN);
    
    % construct right side of C matrix together
    Cright = zeros(2,nL);
    Cright(2,1) = (fso + fsl*ES);
    
    % put together the C matrix
    Param.C(:,:,mc) = [Cleft Cright];
    
    
    % the new 2 state D matrix
    Param.D(:,:,mc) = [(fnl*FN)              0;...
                            0   (fsl*FS)];
end


%% Save the model parameter sets into the base directory

location_for_variables = [base_directory filesep project_name filesep 'MagicC_model_parameter_sets'];
disp(location_for_variables)
save(location_for_variables, 'Param', '-v6'); 

