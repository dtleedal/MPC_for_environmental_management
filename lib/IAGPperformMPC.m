function [Data, emissions] = IAGPperformMPC(year_from_cmd,base_directory,project_name,Data,Param,spin_up_year,save_data_flag,write_SO2_flag)



if isempty(spin_up_year)
    % then the function has not been called with a year of data
    % from the spin up file and therefore the data must be
    % read from in input text file
    % for testing
    %year_from_cmd = 2000,experiment_name = 'exp1';
    [year,...
        this_obs_original,...
        NH_SO2_emission,...
        SH_SO2_emission,...
        f_ghg,...
        f_sw_down,...
        f_sw_up,...
        net_toa,...
        NH_volcano_size,...
        SH_volcano_size]...
        = IAGP_read_input_file(year_from_cmd,base_directory,project_name);
else % instead of reading the text file, unpack the data structure from the spin up
    year = spin_up_year.year;
    this_obs_original = spin_up_year.this_obs_original;
    NH_SO2_emission = spin_up_year.NH_SO2_emission;
    SH_SO2_emission = spin_up_year.SH_SO2_emission;
    f_ghg = spin_up_year.f_ghg;
    f_sw_down = spin_up_year.f_sw_down;
    f_sw_up = spin_up_year.f_sw_up;
    net_toa = spin_up_year.net_toa;
    NH_volcano_size = spin_up_year.NH_volcano_size;
    SH_volcano_size = spin_up_year.SH_volcano_size;
end


if isempty(Param) % model parameters will need to be loaded into the function
    disp('Loading model parameter...')
    try
        load([base_directory filesep project_name filesep 'MagicC_model_parameter_sets'])
    catch
        error('Could not load MagicC model parameters file for project')
    end
    disp('...done')
end
if Param.use_diary_flag == 1
    diary([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'log_file'])
end
disp(['User-specified year: ' num2str(year_from_cmd)])
disp('----')
disp('Loading input file...')


%% get the matrix index from the year
try
    time_index = find(Param.t==year);
catch
    error('The year is outside the simulation period or has not been entered correctly')
end

%% put inputs into timeseries less model-specific baselines
this_obs = zeros(1,3);
this_obs(1,1) = this_obs_original(1) - Param.NH_baseline;     % NH temp perturbation
this_obs(1,2) = this_obs_original(2) - Param.SH_baseline;     % SH temp perturbation
this_obs(1,3) = this_obs_original(3) - Param.sea_ice_baseline;% sea ice extent perturbation



choose_directory = {'state_data' 'backup'}; % the two data directory names
if isempty(Data) % if Data has not been passed into the function
    %% load the data
    disp('Loading data files...')
    
    i = 1; % choose 'state_data'
    load([base_directory filesep Param.project_name filesep choose_directory{i} filesep 'all_data_series'])
    disp('...done')
end % end of should data be loaded

if save_data_flag == 1
    %% back up the data
    
    if Param.do_backup
        disp('Saving data files to backup...')
        i = 2; % choose 'backup'
        save([base_directory filesep Param.project_name filesep choose_directory{i} filesep 'all_data_series'],'Data')
        
        disp('...done')
    end % end of should back up be done
end % end of should Data be saved
%% append other input values
% save the NH and SH temp and sea ice extent
Data.the_observations(time_index,:) = this_obs;
% save the NH and SH volcano size
Data.the_volcanoes(time_index,:) = [NH_volcano_size SH_volcano_size];
% the applied SO2 values
Data.all_NH_SO2(time_index) = NH_SO2_emission;
Data.all_SH_SO2(time_index) = SH_SO2_emission;
%% initially set this function's output to zeros
emissions = zeros(Param.n_model,1);

%% build the fitted_DI_des timeseries if required
if ((Param.t(time_index) == Param.control_start) && (Param.use_fitted_DI_des_flag))
    disp('Creating a fitted DI_des series...')
    % the user wants to make an ice target series
    Data.DI_des = IAGP_build_fitted_DI_des(Data.the_observations(:,3),...
        Param.t,...
        Param.control_start,...
        Param.ice_stabilisation_level,...
        Param.ice_stabilisation_year,...
        Param.sea_ice_baseline);
    
    disp('...done')
end % end build fitted_DIdes



%% calculate integral of error starting the year control is initiated
if (Param.t(time_index) >= Param.control_start)
    Data.all_DE(time_index) = Data.all_DE(time_index-1) + (Data.DI_des(time_index)-this_obs(end));
    
    % apply an integral of error wind up limit
    if Data.all_DE(time_index) < -Param.integral_wind_up_limit
        Data.all_DE(time_index) = -Param.integral_wind_up_limit;
        disp('Integral of error has reached lower limit')
    end
    if Data.all_DE(time_index) > Param.integral_wind_up_limit
        Data.all_DE(time_index) = Param.integral_wind_up_limit;
        disp('Integral of error has reached upper limit')
    end
    
    disp(['Integral of error: ' num2str(Data.all_DE(time_index))])
end

%% big branch statement dependent on whether this is PID_only or not
if ((Param.t(time_index) >= Param.control_start) && (Param.PID_only)) % then we're only doing a simple PI(D) control
    % calculate P and D error series
    Data.all_P_error(time_index) = Data.DI_des(time_index)-this_obs(end);
    Data.all_D_error(time_index) = Data.all_P_error(time_index) - Data.all_P_error(time_index-1);
    RMSD_weighted_SO2 = Param.P_weight.* Data.all_P_error(time_index)...
        + Param.I_weight.* Data.all_DE(time_index)...
        + Param.D_weight.* Data.all_D_error(time_index);
    disp(['Using PID. The PID emissions are: ' num2str(RMSD_weighted_SO2)]) 
    % the result will now replace the RMSD_weighted_SO2 value
    % calculted by MPC when not doing PID control
else % do MPC
    
    %% define optimization settings
    OPTIONS = optimset('maxfunevals',1e10,'maxiter',1e10);
    
    LB = ones(Param.N,1).*0; % optimization lower bound
    UB = ones(Param.N,1).*100; % upper bound
    X0 = ones(Param.N,1)*2; % optimization initial parameter values
    
    %% set up space for Nash-Suttcliff calculations
    Nash_Sut = zeros(Param.n_model,1);
    Nash_Sut_weighted = 0;
    N_S_weighted_SO2 = 0;
    RMSD_weighted_SO2 = 0;
    how_far_back = 40; % how many past values to include in Param.N-S
    
    %% define extrapolation parameters for plotting and forcing regression
    all_the_past_ice = nan(how_far_back+1,Param.n_model); % block to plot
    how_far_back_forcing = 10; % how much past forcing to use for extrapolation
    
    %% define some string formatting string variables used later for printing
    int_format_Nplus1 = '%8d';
    float_format_Nplus1 = '%8.2f';
    for i_count = 1:Param.N
        int_format_Nplus1 = [int_format_Nplus1 '%8d'];      % used in printout of main optimization results
        float_format_Nplus1 = [float_format_Nplus1 '%8.2f'];% used in printout of main optimization results
    end
    float_format_N = '%8.2f'; % used in printout of extrapolated forcing results
    for i_count = 1:Param.N-1
        float_format_N = [float_format_N '%8.2f'];
    end
    float_format_n_model = '%8.2f'; % used in printout of RMSE results
    for i_count = 1:Param.n_model-1
        float_format_n_model = [float_format_n_model '%8.2f'];
    end
    
    
    
    %% transfer inputs to Data.all_forcings matrix
    %  if not use functional estimate of aerosol forcing
    % append the new forcing value onto the end of the forcing time series
    % Data.all_forcings contains the following data:
    %   column  description
    %   1       GHG forcing
    %   2       downwelling short wave
    %   3       upwelling short wave
    %   4       net TOA radiation flux
    %   5       estimated aerosol (+various) forcing
    %   6       NH geoengineering forcing
    %   7       SH geoengineering forcing
    %   8       NH volcanic forcing
    %   9       SH volcanic forcing
    %   10      total NH forcing
    %   11      total SH forcing
    
    Data.all_forcings(time_index,1) = f_ghg;
    Data.all_forcings(time_index,2) = f_sw_down;
    Data.all_forcings(time_index,3) = f_sw_up;
    Data.all_forcings(time_index,4) = net_toa;
    
    %% calculate the non-ghg anthropogenic forcing
    % if year is before instrumental record or
    % if Param.which_f_aero_to_use (the which method switch) is equal to 1
    % use linear relation to GHG forcing
    if ((Param.t(time_index) <= Param.Faero_instruments_begin+1)||(Param.which_f_aero_to_use == 1))  % +1 incase toa etc are NaN up to instrument begins year
        % The GHG forcing is then converted into an estimate of the ne
        f_aero = Param.forcing_conversion(1).*f_ghg + Param.forcing_conversion(2);
        
        disp(['Year is before ' num2str(Param.Faero_instruments_begin) ' using parametetric aerosol forcing function'])
        disp(['Result is: Faero = ' num2str(f_aero)])
        if Param.use_plotting_flag
            H = figure;
            plot(Param.t(1:time_index),[Data.all_forcings(1:time_index-1,5);f_aero],'k')
            title(['Linear estimate of aerosol forcing (' num2str(year(end)) ')'],'fontsize',14)
            xlabel('Year','fontsize',12)
            ylabel('Forcing (\itW m^{-2}\rm)','fontsize',12)
        end % end of do plot
    elseif ((Param.t(time_index) > Param.Faero_instruments_begin+1)&&(Param.which_f_aero_to_use == 2))
        % call to Faero estimation function here
        % get index of Faero instrument start year
        disp('using instrumental estimate of aerosol forcing')
        %try
        
        Faero_start_index = find(Param.t==Param.Faero_instruments_begin);
        total_SO2 = Data.all_NH_SO2(Faero_start_index:time_index,1,1) + Data.all_SH_SO2(Faero_start_index:time_index,1,1);
        
        [H, f_aero] = IAGP_f_aero_estimate(Param.t(Faero_start_index:time_index),...
            Data.all_forcings(Faero_start_index:time_index,:),...
            Data.the_observations(Faero_start_index:time_index,:),...
            total_SO2,...
            Param.two_x_CO2_eq_T,...
            Data.all_a_g_NH(Faero_start_index-1:time_index-1,:),...
            Data.all_a_g_SH(Faero_start_index-1:time_index-1,:),...
            Param.offset,...
            Param.SW_out_initial,...
            Param.phi_bl,...
            Param.use_plotting_flag);
    else
        f_aero = Param.f_aero_pre_specified;
        disp(['Year is after ' num2str(Param.Faero_instruments_begin) ' using pre-defined non-GHG anthropogenic forcing of ' num2str(f_aero)])
        if Param.use_plotting_flag
            H = figure;
            plot(Param.t(1:time_index),[Data.all_forcings(1:time_index-1,5);f_aero],'k')
            title(['pre-defined estimate of aerosol forcing (' num2str(year(end)) ')'],'fontsize',14)
            xlabel('Year','fontsize',12)
            ylabel('Forcing (\itW m^{-2}\rm)','fontsize',12)
        end % end of do plot
    end
    if Param.use_plotting_flag
        set(H,'paperpositionMode','manual','paperunits','normalized','paperposition',[0.05 0.5 0.9 0.5])
        print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_' num2str(year(end))],'-dpsc')
        
        % we can then append each figure to this file as we go using the '-append'
        % flag unfortunately they each go onto a new page :(
        close(H);
    end
    % place the calculated value in column 5 of the Data.all_forcings matrix
    Data.all_forcings(time_index,5) = f_aero;
    
    %% set Data.PI_bl(time_index) and check if update is required
    if Param.use_PI_updating
        if find(Data.PI_update_schedule(:,1) == year,1) > 1
            disp('Performing an update of Param.PI_bl...')
            which_node = find(Data.PI_update_schedule(:,1)==year);
            end_index = find(Param.t == year);
            f = @(g)IAGP_get_ice_err(g, year, Data.PI_update_schedule, which_node, Param.t(1:end_index), Data.the_observations(1:end_index,1), Data.the_observations(1:end_index,3));
            X = lsqnonlin(f,Data.PI_update_schedule(which_node,2),-10,10);% F_aero seems to get too big
            disp(['New value for Param.PI_bl = ' num2str(X)])
            Data.PI_update_schedule(which_node:end,2) = X;
            %Param.PI_bl = X;
        else
            %Param.PI_bl = Data.PI_update_schedule(end,2);
        end
    else
        %Param.PI_bl = Data.PI_update_schedule(1,2); % set Param.PI_bl to the initial value
    end
    disp(['Using value for Data.PI_bl of: ' num2str(Param.PI_bl)])
    %% put the geoengineering inputs in place
    % NH in column 1 SH in column 2
    Data.all_NH_SO2(time_index,1) = NH_SO2_emission;
    Data.all_SH_SO2(time_index,1) = SH_SO2_emission;
    Data.all_forcings(time_index,6:7) = Param.phi_bl.*[NH_SO2_emission SH_SO2_emission];
    
    %% filter the volcano series en-bloc
    % only do this if Param.use_volcano is set otherwise leave volcanic
    % forcing value at 0
    % NH volcanoes
    if Param.use_volcano == 1
    Data.all_forcings(:,8) = Param.phi_bl.*(filter(1-exp(-1/Param.volcanic_residence_time),[1 -exp(-1/Param.volcanic_residence_time)],Data.the_volcanoes(:,1)));
    % SH volcanoes
    Data.all_forcings(:,9) = Param.phi_bl.*(filter(1-exp(-1/Param.volcanic_residence_time),[1 -exp(-1/Param.volcanic_residence_time)],Data.the_volcanoes(:,2)));
    end
        
    
    %% calculate the total forcing
    % NH total forcing = GHG + f_aero + NH geoengineering + NH volcanic
    Data.all_forcings(time_index,10) =...        % total NH forcing
        Data.all_forcings(time_index,1)...       % GHG forcing
        + Data.all_forcings(time_index,5)...     % aerosol (+various)
        + Data.all_forcings(time_index,6)...     % NH geoengineering
        + Data.all_forcings(time_index,8);       % NH volcanic
    % SH total forcing = GHG + f_aero + SH geoengineering + SH volcanic
    Data.all_forcings(time_index,11) =...        % total SH forcing
        Data.all_forcings(time_index,1)...       % GHG forcing
        + Data.all_forcings(time_index,5)...     % aerosol (+various)
        + Data.all_forcings(time_index,7)...     % SH geoengineering
        + Data.all_forcings(time_index,9);       % SH volcanic
    
    %% forcing extrapolation for MPC
    
    if (Param.t(time_index) >= Param.control_start) % do this once controller starts (includes plot and extrapolation)
        % collect forcings to plot common to NH and SH
        % Data.all_forcings (column 1)
        last_few_GHG_forcings_withNaNs = [Data.all_forcings(time_index-how_far_back_forcing:time_index,1);nan(Param.N,1)];
        % Data.all_forcings (column 5)
        last_few_f_aero_forcings_withNaNs = [Data.all_forcings(time_index-how_far_back_forcing:time_index,5);nan(Param.N,1)];
        
        % now do Northern hemisphere
        disp('Extrapolating NH forcing (ignore IRWSM warning)...')
        % total NH forcing - vocanic forcing
        last_few_forcings_NH = Data.all_forcings(time_index-how_far_back_forcing:time_index,10) -...
            Data.all_forcings(time_index-how_far_back_forcing:time_index,8)-Data.all_forcings(time_index-how_far_back_forcing:time_index,6);
        
        last_few_forcings_NH_withNaNs = [last_few_forcings_NH;nan(Param.N,1)];
        
        
        % the low NVR means this is linear interpolation
        smoothed_forcing_NH = irwsm(last_few_forcings_NH_withNaNs,1,1e-5);
        % add the Param.N future NH volcanic forcings to the extrapolation
        % and store in the Data.all_forcings matrix
        Data.all_forcings(time_index+1:time_index+Param.N,10) = smoothed_forcing_NH(end-(Param.N-1):end) + Data.all_forcings(time_index+1:time_index+Param.N,8);
        
        % now do the Southern hemisphere
        disp('Extrapolating SH forcing (ignore IRWSM warning)...')
        % total NH forcing - vocanic forcing
        last_few_forcings_SH = Data.all_forcings(time_index-how_far_back_forcing:time_index,11) -...
            Data.all_forcings(time_index-how_far_back_forcing:time_index,9)-Data.all_forcings(time_index-how_far_back_forcing:time_index,7);
        
        last_few_forcings_SH_withNaNs = [last_few_forcings_SH;nan(Param.N,1)];
        
        
        % the low NVR means this is linear interpolation
        smoothed_forcing_SH = irwsm(last_few_forcings_SH_withNaNs,1,1e-5);
        % add the Param.N future SH volcanic forcings to the extrapolation
        % and store in the Data.all_forcings matrix
        Data.all_forcings(time_index+1:time_index+Param.N,11) = smoothed_forcing_SH(end-(Param.N-1):end) + Data.all_forcings(time_index+1:time_index+Param.N,9);
        disp(['The extrapolated adjusted forcing over the next ' num2str(Param.N) ' years is:'])
        
        
        disp_string = sprintf(float_format_N, Data.all_forcings(time_index+1:time_index+Param.N,10));
        disp(['Norther hemisphere: ' disp_string])
        disp_string = sprintf(float_format_N, Data.all_forcings(time_index+1:time_index+Param.N,11));
        disp(['Southern hemisphere: ' disp_string])
        
        %% forcing extrapolation plot
        if Param.use_plotting_flag
            H = figure;
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),last_few_forcings_NH_withNaNs,'k.','markersize',8)
            hold on
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),last_few_forcings_SH_withNaNs,'r.','markersize',8)
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),last_few_GHG_forcings_withNaNs,'b.','markersize',8)
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),last_few_f_aero_forcings_withNaNs,'g.','markersize',8)
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),smoothed_forcing_NH,'k')
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),smoothed_forcing_SH,'r')
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),Data.all_forcings(time_index-how_far_back_forcing:time_index+Param.N,10),'k:')
            plot(Param.t(time_index-how_far_back_forcing:time_index+Param.N),Data.all_forcings(time_index-how_far_back_forcing:time_index+Param.N,11),'r:')
            
            
            xlabel('Year','fontsize',12)
            ylabel('Forcing (\itW m^{-2}\rm)','fontsize',12)
            legend('NH no volc.','SH no volc.','GHG forcing','aerosol','NH interp.','SH interp.','NH total','SH total')
            title(['Forcing Param.N-step forecasts (' num2str(year) ')'],'fontsize',14)
            set(H,'paperpositionMode','manual','paperunits','normalized','paperposition',[0.05 0.5 0.9 0.5])
            try
                print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_' num2str(year(end))],'-dpsc', '-append')
            catch
                error('Could not print the forcing figure to file result_figure_<year>')
            end
            close(H)
            disp('...done')
        end % end of do plot
    end
    
    
    
    %% the multi model loop
    %
    %   * * * * * * * * * * * * * * * * * * * *
    %   *           beginning of MC           *
    %   %                 loop                *
    %   * * * * * * * * * * * * * * * * * * * *
    %
    for mc = 1:Param.n_model
        
        
        %% state space simulation
        
        
        %    * * * * * * * * * * * * * * * * * * * * * * * * * *
        %    *                                                 *
        %    *                 **IMPORTANT**                   *
        %    *           THE FORCINGS ALL HAVE TO BE           *
        %    *                 MULTIPLIED BY 2                 *
        %    *           BECAUSE OF THE WAY MAGICC             *
        %    *            DEFINES THE HEMISPHERES              *
        %    *                                                 *
        %    * * * * * * * * * * * * * * * * * * * * * * * * * *
        
        % extract NH and SH forcings and multiply by 2
        NH_SH_forcing_times_2 = Data.all_forcings(time_index,10:11)'.*2;
        %
        % The function template is:
        % [X_k, Y_k] = stateSpaceSim(X_k_1,U_k,Param.dt,Param.A,Param.B,Param.C,Param.D)
        [last_state, Y_hat_k] = stateSpaceSim(Data.all_X_states(:,time_index-1,mc), NH_SH_forcing_times_2, Param.dt,...
            Param.A(:,:,mc),Param.B(:,:,mc),Param.C(:,:,mc),Param.D(:,:,mc));
        
        
        Data.all_X_states(:,time_index,mc) = last_state;
        % place results into matrix
        Data.all_Y_states(:,time_index,mc) = Y_hat_k;
        
        
        
        %% adaptive gain
        
        if Param.use_T_adaptive_gain_flag
            % NH temperature
            g_k_1_NH = Data.all_a_g_NH(time_index-1,mc);
            p_k_1_NH = Data.all_p_k_NH(time_index-1,mc);
            y_hat_NH = Y_hat_k(1);%.*g_k_1_NH; % NH temp * gain
            y = Data.the_observations(time_index,1);%NH temp
            [g_k_NH, p_k_NH] = IAGP_adaptive_gain( g_k_1_NH, p_k_1_NH, y_hat_NH+Param.offset, y+Param.offset, Param.a_g_NVR_T);
            
            Data.all_a_g_NH(time_index,mc) = g_k_NH;
            Data.all_p_k_NH(time_index,mc) = p_k_NH;
            
            
            % SH temperature
            g_k_1_SH = Data.all_a_g_SH(time_index-1,mc);
            p_k_1_SH = Data.all_p_k_SH(time_index-1,mc);
            y_hat_SH = Y_hat_k(2);%.*g_k_1_SH; % SH temp * gain
            y = Data.the_observations(time_index,2); % SH temp
            [g_k_SH, p_k_SH] = IAGP_adaptive_gain( g_k_1_SH, p_k_1_SH, y_hat_SH+Param.offset, y+Param.offset, Param.a_g_NVR_T);
            
            Data.all_a_g_SH(time_index,mc) = g_k_SH;
            Data.all_p_k_SH(time_index,mc) = p_k_SH;
        else
            % set both NH and SH adaptive gains to 1
            Data.all_a_g_NH(time_index,mc) = 1;
            Data.all_a_g_SH(time_index,mc) = 1;
            g_k_NH = 1;
            g_k_SH = 1;
        end
        if Param.use_SI_adaptive_gain_flag
            g_k_1_SI = Data.all_a_g_SI(time_index-1,mc);
            p_k_1_SI = Data.all_p_k_SI(time_index-1,mc);
            y_hat_SI = ((((((Y_hat_k(1)+Param.offset).*g_k_NH)-Param.offset)).*Param.PI_bl)); % NH temp * NH_gain * SI_gain * PI_bl
            y = Data.the_observations(time_index,3);
            [g_k_SI, p_k_SI] = IAGP_adaptive_gain( g_k_1_SI, p_k_1_SI, y_hat_SI+Param.offset, y+Param.offset, Param.a_g_NVR_SI);
            Data.all_a_g_SI(time_index,mc) = g_k_SI;
            Data.all_p_k_SI(time_index,mc) = p_k_SI;
        else
            % set sea ice adaptive gain to 1
            Data.all_a_g_SI(time_index,mc) = 1;
            g_k_SI = 1;
        end
        
        
        
%         %% fill out the model error matrix (not used at the moment)
%         % SE for NH and SH temperature
%         Data.model_errors(time_index,1:2,mc) = (this_obs(1:2) - Data.all_Y_states(1:2,time_index,mc)').^2;
%         % time-weighted SSE for NH and SH
%         Data.model_errors(time_index,3:4,mc) = Param.SSE_weight.*Data.model_errors(time_index-1,3:4,mc)...
%             + Data.model_errors(time_index,1:2,mc);
        
        %% Adjust all_Y_states by adaptive gain amount
        %     Data.all_Y_states(1,time_index,mc) = g_k_NH.*Data.all_Y_states(1,time_index,mc);
        %     Data.all_Y_states(2,time_index,mc) = g_k_SH.*Data.all_Y_states(2,time_index,mc);
        
        %% convert MagicC output to sea ice extent and store
        Data.all_ice_states(time_index, 1, mc) = (((((Data.all_Y_states(1,time_index,mc)+Param.offset).*g_k_NH)-Param.offset).*Param.PI_bl)+Param.offset.*g_k_SI)-Param.offset;
        
        %% MPC optimization
        if (Param.t(time_index) >= Param.control_start)
            disp('-------------------------')
            disp(['Summary for model ' num2str(mc)])
            disp(['Target value:     ' num2str(Data.DI_des(time_index))])
            % at this point Param.PI_bl is not re-estimated
            disp(['Observed estimate of Sea ice: ' num2str(Data.all_ice_states(time_index, 1, mc))])
            
            Xpars = fminsearchbnd('IAGP_optimize_SO2',X0,LB,UB,OPTIONS,...
                Param.A(:,:,mc),Param.B(:,:,mc),Param.C(:,:,mc),Param.D(:,:,mc),...
                Data.all_forcings(:,10:11).*2,... % NH and SH total forcing
                Data.all_NH_SO2.*2,...
                Data.all_X_states(:,:,mc),...
                Data.all_Y_states(:,:,mc),...
                Param.optimize_weights,...
                Data.DI_des,... % the target
                Data.all_DE,... % the integral of error
                g_k_NH,...
                g_k_SH,...
                g_k_SI,...
                Param.offset,...
                Param.phi_bl,...
                Param.PI_bl,...
                Param.use_last_observation,...
                Param.N,time_index,Param.dt);
            %catch
            % error('Failured during call to optimization')
            %end
            % spin out the forecasts
            % ----------------------
            %try
            the_ice_single = IAGP_get_n_step_ice...
                (Xpars,...% these are the SO2 emissions (already *2)
                Param.A(:,:,mc),Param.B(:,:,mc),Param.C(:,:,mc),Param.D(:,:,mc),...
                Data.all_forcings(:,10:11).*2,... % NH and SH total forcing
                Data.all_X_states(:,:,mc),...
                Data.all_Y_states(:,:,mc),...
                g_k_NH,...
                g_k_SH,...
                g_k_SI,...
                Param.offset,...
                Param.phi_bl,...
                Param.PI_bl,...
                Param.N,time_index,Param.dt);
            Data.all_ice_states(time_index,2:end,mc) = the_ice_single;
            %catch
            %    error('Failed during call to Data.all_ice_states')
            %end
            % ----------------------
            half_Xpars = Xpars./2; % all optimal emissions are twice the real size
            % due to MagicC splitting the forcing
            
            %% make text output formatting strings
            
            
            % output information to the screen
            text_string_year = sprintf(['Year:               ' int_format_Nplus1], Param.t(time_index:time_index+Param.N));
            disp(text_string_year)
            
            text_string0 = sprintf(['Assumed forcing:    ' float_format_Nplus1], Data.all_forcings(time_index:time_index+Param.N,10));
            disp(text_string0)
            text_string1 = sprintf(['Optimal emissions:  ' float_format_Nplus1], [Data.all_NH_SO2(time_index,1) half_Xpars']);
            disp(text_string1)
            text_string2 = sprintf(['Modelled ice path:  ' float_format_Nplus1], [this_obs(end) the_ice_single]);
            disp(text_string2)
            text_string3 = sprintf(['Target ice path:    ' float_format_Nplus1], Data.DI_des(time_index:time_index+Param.N));
            disp(text_string3)
            error_path = [Data.all_DE(time_index) (Data.DI_des(time_index+1:time_index+Param.N)'-the_ice_single)];
            text_string4 = sprintf(['Modelled error:     ' float_format_Nplus1], [(Data.DI_des(time_index) - this_obs(end)) error_path(2:end)]);
            disp(text_string4)
            text_string5 = sprintf(['Modelled IofE:      ' float_format_Nplus1], cumsum(error_path));
            disp(text_string5)
            Data.all_NH_SO2(time_index+1,2:end,mc) = half_Xpars;% this is the emissions for the real-word i.e., MagicC-derived divided by 2
            
            
            
            
            
        end % end of "are we in control period?" conditional
        
        
        
    end % end of mc ensemble loop
    
%     %% fill out relative SSE score for each model (not used at the moment
%     % NH
%     term_1 = 1 - Data.model_errors(time_index,3,:)./sum(Data.model_errors(time_index,3,:));
%     NH_model_errors = term_1./sum(term_1);
%     % SH
%     term_1 = 1 - Data.model_errors(time_index,4,:)./sum(Data.model_errors(time_index,4,:));
%     SH_model_errors = term_1./sum(term_1);
%     
%     Data.model_errors(time_index,5,:) = mean([NH_model_errors SH_model_errors]);
    %% plotting section
    if Param.use_plotting_flag
        if (Param.t(time_index) >= Param.control_start)
            H = figure;
            % NH and SH temperature plot
            plot(Param.t(time_index-how_far_back:time_index),Data.the_observations(time_index-how_far_back:time_index,1:2),'linewidth',2)
            hold on
            % estimate (with data assimilation)
            for mc = 1:Param.n_model
                plot(Param.t(time_index-how_far_back:time_index),Data.all_Y_states(1:2,time_index-how_far_back:time_index,mc))
            end
            xlabel('Year')
            ylabel('Temperature (perturbation)')
            title({'Observed NH temperature (thick blue), model ensemble NH temperature (blue)'; 'Observed SH temperature (thick green), model ensemble SH temperature (green)'})
            try
                print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_' num2str(year(end))],'-dpsc', '-append')
            catch
                error('Could not print the temperature figure to file NH_SH_temp_<year>')
            end
            close(H)
            
            % emissions plot
            H = figure;
            plot(Param.t(find(Param.t==(Param.control_start-5)):time_index),Data.all_NH_SO2(find(Param.t==(Param.control_start-5)):time_index,1,1),'linewidth',2)
            hold on
            for mc = 1:Param.n_model
                plot(Param.t(time_index+1:time_index+Param.N),Data.all_NH_SO2(time_index,2:end,mc),'linewidth',1)
            end
            xlabel('Year')
            ylabel('NH SO2 emissions (annual)')
            title('NH SO2 emissions + forecasts')
            try
                print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_' num2str(year(end))],'-dpsc', '-append')
            catch
                error('Could not print the emissions figure to file result_figure_<year>')
            end
            close(H)
            
            % ice plot
            % make this start at 2000
            index_2000 = find(Param.t==2000);
            H = figure;
            plot(Param.t(index_2000:time_index),Data.the_observations(index_2000:time_index,3),'linewidth',2)
            hold on
            for mc = 1:Param.n_model
                plot(Param.t(index_2000:time_index),Data.all_ice_states(index_2000:time_index,1,mc))
                plot(Param.t(time_index+1:time_index+Param.N),Data.all_ice_states(time_index,2:end,mc),'r')
            end
            
            % add the target line
            plot(Param.t(index_2000:time_index+1),Data.DI_des(index_2000:time_index+1),'k','linewidth',2)
            xlabel('Year')
            ylabel('Min sea ice extent (perturbation)')
            title({'Observed sea ice minimum (thick blue), target (thick black),'; 'model ensemble 1-step forecasts (thin blue), Param.N-step forecast (red)'})
            try
                print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_' num2str(year(end))],'-dpsc', '-append')
            catch
                error('Could not print the sea ice figure to file result_figure__<year>')
            end
            close(H)
            
        end
    end % end of do plot
    
    
    %% RMSE calculation section
    %
    disp('Calculating RMSD...')
    if (Param.t(time_index) >= Param.control_start)
        RMSD_model_weights = ones(Param.n_model,1);
        if ((Param.use_weighting) && (Param.n_model > 1))
            
            RMSD_value = ones(Param.n_model,1);
            if (Param.t(time_index) > Param.control_start)
                % observation_set is all the sea ice observations (less baseline)
                % from one year after the start of the control period to the present time step
                observation_set = Data.the_observations(find(Param.t==Param.control_start+1):time_index,end);
                exponential_weights = flipud(filter(1,[1 -Param.RMSD_exponential_weight],[1; zeros(length(observation_set)-1,1)]));
                for RMSE_model_index = 1:Param.n_model
                    % one_step_forecast_set is the one year ahead forecast made by one of the 10
                    % models for the sea ice (less baseline) made from the start of the control period
                    % to one year before the present time step
                    one_step_forecast_set = Data.all_ice_states(find(Param.t==Param.control_start):time_index-1,2,RMSE_model_index);
                    weighted_sea_ice_error_set = (observation_set - one_step_forecast_set).*exponential_weights;
                    RMSD_value(RMSE_model_index) = sqrt(sum(weighted_sea_ice_error_set.^2)./length(observation_set));
                end
                RMSD_model_weights = 1 - (RMSD_value./max(RMSD_value));
                RMSD_model_weights = RMSD_model_weights./max(RMSD_model_weights);
                RMSD_model_weights = RMSD_model_weights.^3;
                
                
            end
        end
        
        sum_of_RMSD_model_weights = sum(RMSD_model_weights);
        RMSD_model_weights = RMSD_model_weights./sum_of_RMSD_model_weights;
        weighted_sum_of_SO2 = sum(RMSD_model_weights.*reshape(Data.all_NH_SO2(time_index+1,2,:),Param.n_model,1)); % Param.n_model is the number of MagicC models
        sum_of_RMSD_weights = sum(RMSD_model_weights);
        RMSD_weighted_SO2 = weighted_sum_of_SO2./sum_of_RMSD_weights;
        
        text_string = sprintf(['RMSD:              ' float_format_n_model], RMSD_model_weights);
        disp(text_string)
        text_string = sprintf(['each SO2:          ' float_format_n_model], reshape(Data.all_NH_SO2(time_index+1,2,:),Param.n_model,1));
        disp(text_string)
        text_string = sprintf('RMSD weighted SO2: %8.2f', RMSD_weighted_SO2);
        disp(text_string)
    end
    
end % end of do PID_only conditional branch
%% plant failure section
disp('  ')
disp('----------------------------------------------')
disp(['                ' num2str(Param.t(time_index)) ' REPORT                '])
disp('----------------------------------------------')



if (Param.t(time_index) >= Param.control_start) % check for plant failure
    if Param.use_plant_failure == 1
        % plant failure module
        % --------------------
        % check Data.the_failures time series to see how long it has been
        % since the last failure
        % Data.the_failures is a vector the same length as Param.t with a value
        % of 0 (no failure) 1 (January failure) 2 (february failure) or
        % 3 (march failure)
        when_failed = find(Data.the_failures~=0);
        if ~isempty(when_failed) % there has been a failure since 2019
            timeSinceFix = Param.t(time_index) - Param.t(when_failed(end));
        else % there hasn'Param.t been a failure since 2019
            timeSinceFix = Param.t(time_index) - Param.control_start;
        end
        
        disp('  ')
        disp(['we are at the end of year ' num2str(Param.t(time_index)) '. It''s been ' num2str(timeSinceFix) ' years since a plant failure'])
        
        probOfSurvival = 1 - wblcdf(timeSinceFix,Param.weibul_a,Param.weibul_b);
        
        flip = rand;
        if flip > probOfSurvival
            which_month = ceil(rand*3);
            switch which_month
                case 1
                    disp('Param.A PLANT FAILURE WILL PUT FIRST 7 WEEKS IN WITH WEEKS 8 TO 14!')
                    Data.the_failures(time_index) = 1;
                case 2
                    disp('Param.A PLANT FAILURE WILL PUT WEEKS 8 TO 14 IN WITH WEEKS 15 TO 21!')
                    Data.the_failures(time_index) = 2;
                case 3
                    disp('Param.A PLANT FAILURE WILL PUT WEEKS 15 TO 21 IN WITH WEEKS 8 TO 14!')
                    Data.the_failures(time_index) = 3;
            end
            
            
        else
            disp('The plant will not fail next year')
            Data.the_failures(time_index) = 0;
        end
    else
        disp('Not using plant failure mode')
        Data.the_failures(time_index) = 0;
    end
    
end





% new weekly emissions breakdown
% Week    SO2emissionsmultiplier
%     1.0000    0.0257
%     2.0000    0.0328
%     3.0000    0.0396
%     4.0000    0.0468
%     5.0000    0.0537
%     6.0000    0.0593
%     7.0000    0.0647
%     8.0000    0.0680
%     9.0000    0.0709
%    10.0000    0.0711
%    11.0000    0.0709
%    12.0000    0.0683
%    13.0000    0.0642
%    14.0000    0.0592
%    15.0000    0.0533
%    16.0000    0.0461
%    17.0000    0.0381
%    18.0000    0.0297
%    19.0000    0.0209
%    20.0000    0.0124
%    21.0000    0.0041


% break down emissions into 10 blocks

% the 21 weights determined by AJ
SO2_weights = [0.0257 0.0328 0.0396 0.0468 0.0537 0.0593 0.0647 0.0680 ...
    0.0709 0.0711 0.0709 0.0683 0.0642 0.0592 0.0533 0.0461 0.0381 ...
    0.0297 0.0209 0.0124 0.0041];
% the sum of weights of each section
mean_section_weights(1) = mean(SO2_weights(1:7));
mean_section_weights(2) = mean(SO2_weights(8:14));
mean_section_weights(3) = mean(SO2_weights(15:21));

switch Data.the_failures(time_index)
    case 0 % no failure
        section_1_emissions = RMSD_weighted_SO2.*SO2_weights(1:7);
        section_2_emissions = RMSD_weighted_SO2.*SO2_weights(8:14);
        section_3_emissions = RMSD_weighted_SO2.*SO2_weights(15:21);
        
    case 1 % section 1 failure
        section_1_emissions = 0.*SO2_weights(1:7);
        section_2_emissions = RMSD_weighted_SO2.*(SO2_weights(8:14)+mean_section_weights(1));
        section_3_emissions = RMSD_weighted_SO2.*SO2_weights(15:21);
        
    case 2 % section 2 failure
        section_1_emissions = RMSD_weighted_SO2.*SO2_weights(1:7);
        section_2_emissions = 0.*SO2_weights(8:14);
        section_3_emissions = RMSD_weighted_SO2.*(SO2_weights(15:21)+mean_section_weights(2));
        
    case 3 % section 3 failure
        section_1_emissions = RMSD_weighted_SO2.*SO2_weights(1:7);
        section_2_emissions = RMSD_weighted_SO2.*(SO2_weights(8:14)+mean_section_weights(3));
        section_3_emissions = 0.*SO2_weights(15:21);
    otherwise
        disp('*** ERROR *** THE FAILURES SWITCH IS NOT WORKING!!!')
end
all_SO2_emissions = [section_1_emissions section_2_emissions section_3_emissions];

if write_SO2_flag % during spin up processing there is no point writing the SO2 file
    % make SO2 output file whith SO2 emissions to use for each week
    IAGP_make_SO2_output_file(...
        base_directory,...
        Param.project_name,...
        year+1,...
        all_SO2_emissions)
end
if save_data_flag == 1
    disp('Saving data files...')
    i = 1; % i = 1 sets choose_directory equal to the string 'state_data'
    save([base_directory filesep Param.project_name filesep choose_directory{i} filesep 'all_data_series'],'Data')
    
    disp('...done')
    
end

emissions = RMSD_weighted_SO2;
% switch the diary keeping off (if it has been specified to be on by the user)
if Param.use_diary_flag == 1
    diary OFF
end



