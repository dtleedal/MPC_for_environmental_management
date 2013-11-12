function Data = IAGP_data_setup_ONE_OFF(base_directory,project_name)
disp(' ')
disp(' ')
disp('* * * * * * * * * * * * * * * * * * * * * * * *')
disp('*                                             *')
disp('*                W A R N I N G                *')
disp('*          this will overwrite existing       *')
disp('*         results and should only be used     *')
disp('*        once to initialize the experiment    *')
disp('*                                             *')
reply = input('*        DO YOU WANT TO CONTINUE (Y/N):       *\n*                                             *\n* * * * * * * * * * * * * * * * * * * * * * * *\n>','s');

if (reply~='Y' && reply~='y')
    disp('Leaving the function without writing any data.')
    return
else
    disp(' ')
    disp('OK then. Writing initialization files to:')
    disp([base_directory filesep project_name filesep 'state_data'])
end


%%

% this script sets up the data structures that will store the walk through
% results and sets all values to zero.
%
% * * * * * * * * * * * * * * * * * * * * * * * *
% *                                             *
% *                W A R N I N G                *
% *              do not run during              *
% *          the experiment as it will          *
% *        overwrite the existing results!      *
% *                                             *
% * * * * * * * * * * * * * * * * * * * * * * * *

% get the model variables
% define the base working directory this will end with the directory: \MPC_for_Lawrence
%base_directory = 'C:\Documents and Settings\leedald\My Documents\Dropbox\IAGP_HadGEM_control\AJs_control_model\MPC_for_Lawrence';
%base_directory = '/Users/davidleedal/Dropbox/IAGP_HadGEM_control/AJs_control_model/MPC_for_Lawrence';

% retrieve the model parameters used to size the result arrays
location_of_variables = [base_directory filesep project_name filesep 'MagicC_model_parameter_sets'];
load(location_of_variables)


% Make space for and save variables to hold
%

Data.all_X_states = zeros(2*Param.nL,length(Param.t),Param.n_model);
Data.all_Y_states = zeros(2,length(Param.t),Param.n_model);
Data.all_ice_states = nan(length(Param.t),Param.N+1,Param.n_model); % store sea ice estimate in col 1 + N forecasts
Data.all_p_k_NH = zeros(length(Param.t),Param.n_model); % the error variance for the NH temperature
Data.all_p_k_NH(1,:) = Param.P_k_0;
Data.all_p_k_SH = zeros(length(Param.t),Param.n_model); % the error variance estimate for the SH temperature
Data.all_p_k_SH(1,:) = Param.P_k_0;
Data.all_p_k_SI = zeros(length(Param.t),Param.n_model); % the error variance estimate for the adaptive gains
Data.all_p_k_SI(1,:) = Param.P_k_0;
Data.all_a_g_NH = zeros(length(Param.t),Param.n_model); % the adaptive gain for the NH temperature
Data.all_a_g_NH(1,:) = 0;
Data.all_a_g_SH = zeros(length(Param.t),Param.n_model); % the adaptive gain for the SH temperature
Data.all_a_g_SH(1,:) = 0;
Data.all_a_g_SI = zeros(length(Param.t),Param.n_model); % the adaptive gain for the sea ice
Data.all_a_g_SI(1,:) = 1;
Data.all_NH_SO2 = zeros(length(Param.t),Param.N+1);% NH geoengineering emissions (col 1) + N forecasts
Data.all_SH_SO2 = zeros(length(Param.t),Param.N+1);% SH geoengineering emissions (col 1) + N forecasts
Data.all_DE = zeros(length(Param.t),1); % integral of error
if Param.PID_only
    Data.all_P_error = zeros(length(Param.t),1); %  error
    Data.all_D_error = zeros(length(Param.t),1); %  derivative of error
end
Data.the_volcanoes = zeros(length(Param.t),2); % NH and SH volcanic emissions series
Data.all_forcings = zeros(length(Param.t),11);
Data.the_failures =  zeros(length(Param.t),1);
Data.the_observations =  zeros(length(Param.t),3);
if Param.use_fitted_DI_des_flag ~= 1
    load([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'sea_ice_target'])
    Data.DI_des = [nan; sea_ice_target] - Param.sea_ice_baseline;
else
    Data.DI_des = nan(length(Param.t),1);
end
% make the PI parameter series
PI_update_schedule_year = [Param.start_year; (Param.control_start:10:Param.end_year)'];
PI_update_N_values = ones(length(PI_update_schedule_year),1).*Param.PI_bl;
Data.PI_update_schedule = [PI_update_schedule_year PI_update_N_values];
Data.model_errors = zeros(length(Param.t),6,Param.n_model); % (1) NH sse, (2) SH sse, (3) NH weighted SSE,
Data.model_errors(1,5:7,:) = 1/Param.n_model;               % (4) SH weighted SSE, (5) NH relative score,                                                        % error, NH
% (6) SH relative score

choose_directory = {'state_data' 'backup'};
for i = 1:2
    save([base_directory filesep project_name filesep choose_directory{i} filesep 'all_data_series'],'Data')
    
    % * * * * * * * * * * * * * * * * * * * * * * * *
    % *                                             *
    % *              I M P O R T A N T              *
    % *            the same baseline that is        *
    % *      subtracted from the sea ice series     *
    % *       must be subtracted from the target!   *
    % *                     5.5                     *
    % * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
end


