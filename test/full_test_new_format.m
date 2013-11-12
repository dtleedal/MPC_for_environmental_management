fclose('all')
clear all
close all
rng(55);

%[A,B,C,D,L] = make_test_model;
load The_test_model_for_paper
[L, tcr] = IAGP_get_lambda_and_tcr( A,B,C,D,40,0.1 )
% matrices for test model results
% parameters for specific model
nL = 40;
t = (1859:2110)';
all_X_states = zeros(2*nL,length(t));
all_Y_states = zeros(2,length(t));
all_sea_ice = zeros(length(t),1);
albedo_test = 0.29;
PI_bl_test = -2.1;
phi_bl_test = -0.16;
% get RCP4.5 GHG forcing 1859 to 2110
GHG_forcing = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','F155:F406');
% total_anth = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','E155:E406');
% total_all = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','B155:B406');
% volc_forcing = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','C155:C406');
% volc_emis = (1./Param.phi_bl).*(filter([1 -exp(-1/Param.volcanic_residence_time)],1-exp(-1/Param.volcanic_residence_time),volc_forcing));
% %dir_aero = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','AP155:AP406');
% %volc_forcing = zeros(size(volc_forcing));
% delta_SI = xlsread([pwd filesep 'test' filesep 'RCP45_MIDYEAR_RADFORCING.xlsx'],'RCP45_MIDYEAR_RADFORCING','D155:D406');
% dir_aero = total_anth - GHG_forcing
SI_in_mean = 341.5; % average sounding TOA input from sun
%test_so2 = zeros(252,1);test_so2(152:end) = 5;

% make a set of NH and SH volcanic events
volcano_series = zeros(length(t),2);
volcano_series(170,1) = 5; % 2028 NH hemisphere
volcano_series(177,2) = 9; % 2035 SH hemisphere
volcano_series(190,1) = 8; % 2048 NH hemisphere
volcano_series(192,2) = 7; % 2050 SH hemisphere
volcano_series(202,1) = 5; % 2060 NH hemisphere
% make a synthetic F_aero pattern
F_aero = zeros(length(t),1);
F_int1 = -0.4.*GHG_forcing(1:150)-0.115;
F_int2 = linspace(F_int1(end),0,102)';
F_int3 = F_int2 + (randn(size(F_int2)).*0);
F_int4 = irwsm([F_int1; F_int3] ,0,1e-2);
F_aero = F_int4;
% plot(t,F_aero)
% make a synthetic Incoming SW pattern
F_int1 = randn(length(t),1).*0;
SW_in = irwsm(F_int1 ,0,1e-1)+SI_in_mean;
SW_in(1:2) = SI_in_mean;
%plot(t,SW_in)
% make a synthetic albedo pattern
F_int1 = randn(length(t),1).*0;
albedo = irwsm(F_int1 ,0,1e-2)+albedo_test;
albedo(1:2) = albedo_test; % initial unperturbed conditions
%plot(albedo)

% use SW_in and albedo to produce SW_out
SW_out = SW_in.*albedo;
%plot(SW_out)
noise_1_b = randn*0.08;
noise_2_b = randn*0.08;
frac_of_b = 0.3;
NH_v_f = 0;
SH_v_f = 0;
volcanic_residence_time = 1;
NH_baseline = 13.1115;
SH_baseline = 13.5478;
global_baseline = 13.3297;
sea_ice_baseline = 5.5000;
Faero_instruments_begin = 1990;
dt = 0.1;
spin_up_data = zeros((2017-1)-1860,12);
row_index = 1;
for the_year = t(2):2016 % year before control begins
    
    time_ind = find(t == the_year)
    year = t(time_ind);
    %for the_year = 1991:end_year
    
    Fghg = GHG_forcing(time_ind);
    
    
    % get SO2 emissions
    NH_SO2 = 0;
    Fso2 = 0;
    % enough to run the model
    
    % make a chance of a random NH or SH volcano
    %     if rand < 0.08; volcano_size = rand.*12;else volcano_size = 0;end
    %     if rand<0.5; NH_v = volcano_size;SH_v = 0;else SH_v = volcano_size;NH_v = 0;end
    NH_v = volcano_series(time_ind,1);
    SH_v = volcano_series(time_ind,2);
    
    NH_v_f = (1-exp(-1/volcanic_residence_time)).*NH_v + exp(-1/volcanic_residence_time).* NH_v_f;
    SH_v_f = (1-exp(-1/volcanic_residence_time)).*SH_v + exp(-1/volcanic_residence_time).* SH_v_f;
    
    %Data.all_forcings(:,8) = Param.phi_bl.*(filter(1-exp(-1/Param.volcanic_residence_time),[1 -exp(-1/Param.volcanic_residence_time)],Data.the_volcanoes(:,1)));
    TOA_NH = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + Fso2 ...
        + phi_bl_test.*NH_v_f ...
        - (SW_out(time_ind) - SW_out(1));
    
    TOA_SH = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + phi_bl_test.*SH_v_f ...
        - (SW_out(time_ind) - SW_out(1));
    U_k = [TOA_NH; TOA_SH].*2;
    [X_k, Y_k] = stateSpaceSim(all_X_states(:,time_ind-1),U_k,dt,A,B,C,D);
    %     T_noise = randn*0;
    %     Y_k = Y_k + T_noise;
    all_X_states(:,time_ind) = X_k;
    all_Y_states(:,time_ind) = Y_k;
    FT = (3.7./L).*mean(all_Y_states(:,time_ind));
    TOA = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + Fso2./2 ...
        + phi_bl_test.*NH_v_f ...
        + phi_bl_test.*SH_v_f ...
        -FT ...
        - (SW_out(time_ind) - SW_out(1));
    NH_t = Y_k(1) + NH_baseline;
    SH_t = Y_k(2) + SH_baseline;
    all_sea_ice(time_ind) = (PI_bl_test .* Y_k(1))+sea_ice_baseline;
    sea_i = all_sea_ice(time_ind);
    
    
    if year > Faero_instruments_begin
        f_sw_down = SW_in(time_ind);
        f_sw_up = SW_out(time_ind);
        
    else
        f_sw_down = NaN;
        f_sw_up = NaN;
        TOA = NaN;
    end
    
    
    SH_SO2 = 0;
    % sum SO2 emissions from IAGP file
    % calculate sea ice min using PI_b
    noise_1_a = randn*0.08;noise_2_a = randn*0.08;
    
    spin_up_data(row_index,:) = [year ...
        NH_t+(frac_of_b*noise_1_b+(1-frac_of_b)*noise_1_a) ...
        SH_t+(frac_of_b*noise_2_b+(1-frac_of_b)*noise_2_a) ...
        sea_i+(randn*0.2) ...
        NH_SO2 ...
        SH_SO2 ...
        Fghg ...
        f_sw_down ...
        f_sw_up ...
        TOA ...
        NH_v ...
        SH_v];
    row_index = row_index + 1;
end
save spin_up_data spin_up_data


% testing script for MPC_for_Ben format
% make a test model
% set up directories variables etc
project_name = 'single_model_no_volcanic'
IAGP_directory_setup(pwd,project_name)
copyfile([pwd filesep 'spin_up_data.mat'],...
    [pwd filesep, project_name filesep 'inputs_and_outputs'])
Param = IAGP_model_setup(pwd,project_name)
% copy sea ice target to new project
% copyfile([pwd filesep 'sea_ice_targets' filesep 'democratic_sea_ice_target.mat'],...
%     [pwd filesep, project_name filesep 'inputs_and_outputs'])
Data = IAGP_data_setup_ONE_OFF(pwd,project_name)
[Data, Param] = IAGP_batch_process_spin_up_data(pwd, project_name, Data, Param)



write_SO2_flag = 1;
save_data_flag = 1;
for the_year = 2017:2100
    %for the_year = 2018:Param.end_year
    time_ind = find(t == the_year);
    year = t(time_ind);
    %for the_year = 1991:end_year
    
    Fghg = GHG_forcing(time_ind);
    
    
    % get SO2 emissions
    NH_SO2 = read_SO2_file(pwd, project_name, year);
    Fso2 = phi_bl_test.*NH_SO2;
    % enough to run the model
    
    % make a chance of a random NH or SH volcano
    %     if rand < 0.08; volcano_size = rand.*12;else volcano_size = 0;end
    %     if rand<0.5; NH_v = volcano_size;SH_v = 0;else SH_v = volcano_size;NH_v = 0;end
    NH_v = volcano_series(time_ind,1);
    SH_v = volcano_series(time_ind,2);
    
    NH_v_f = (1-exp(-1/volcanic_residence_time)).*NH_v + exp(-1/volcanic_residence_time).* NH_v_f;
    SH_v_f = (1-exp(-1/volcanic_residence_time)).*SH_v + exp(-1/volcanic_residence_time).* SH_v_f;
    
    %Data.all_forcings(:,8) = Param.phi_bl.*(filter(1-exp(-1/Param.volcanic_residence_time),[1 -exp(-1/Param.volcanic_residence_time)],Data.the_volcanoes(:,1)));
    TOA_NH = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + Fso2 ...
        + phi_bl_test.*NH_v_f ...
        - (SW_out(time_ind) - SW_out(1));
    
    TOA_SH = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + phi_bl_test.*SH_v_f ...
        - (SW_out(time_ind) - SW_out(1));
    U_k = [TOA_NH; TOA_SH].*2;
    [X_k, Y_k] = stateSpaceSim(all_X_states(:,time_ind-1),U_k,dt,A,B,C,D);
    %     T_noise = randn*0;
    %     Y_k = Y_k + T_noise;
    all_X_states(:,time_ind) = X_k;
    all_Y_states(:,time_ind) = Y_k;
    FT = (3.7./L).*mean(all_Y_states(:,time_ind));
    TOA = GHG_forcing(time_ind) ...
        + F_aero(time_ind) ...
        + Fso2./2 ...
        + phi_bl_test.*NH_v_f ...
        + phi_bl_test.*SH_v_f ...
        -FT ...
        - (SW_out(time_ind) - SW_out(1));
    NH_t = Y_k(1) + NH_baseline;
    SH_t = Y_k(2) + SH_baseline;
    all_sea_ice(time_ind) = (PI_bl_test .* Y_k(1))+sea_ice_baseline;
    sea_i = all_sea_ice(time_ind);
    
    
    
    f_sw_down = SW_in(time_ind);
    f_sw_up = SW_out(time_ind);
    
    
    
    SH_SO2 = 0;
    % sum SO2 emissions from IAGP file
    % calculate sea ice min using PI_b
    noise_1_a = randn*0.08;noise_2_a = randn*0.08;
    
    test_write_input_file( pwd, project_name, year, NH_t+(frac_of_b*noise_1_b+(1-frac_of_b)*noise_1_a), SH_t+(frac_of_b*noise_2_b+(1-frac_of_b)*noise_2_a), sea_i+(randn*0.2), NH_SO2, SH_SO2,...
        Fghg, f_sw_down, f_sw_up, TOA, NH_v, SH_v)
    noise_1_b = noise_1_a;noise_2_b = noise_2_a;
    [Data, emissions] = IAGPperformMPC(the_year, pwd, project_name, [], [], [], save_data_flag, write_SO2_flag);
    
end

fclose('all')





% continue from here with
IAGP_eazy_plot(pwd, Data, Param, 1980, 2100, F_aero)


% results plots for paper
which_project = {'simple multi-model inversion', 'multi-model with integral of error', 'adaptive multi-model with integral of error', 'adaptive multi-model inversion','Proportional, Integral Derivative (PID)' 'single_model_no_volcanic'};
plot_save_dir = 'C:\Documents and Settings\leedald\My Documents\Dropbox\Sea_ice_control_papers\policy_as_recursive_process\PRP_figures'
plot_save_dir = '/Users/davidleedal/Dropbox/Sea_ice_control_papers/policy_as_recursive_process/PRP_figures'

for project_no = 1:5
    
    
    % load each projects data
    switch project_no
        case 1
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 2
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 3
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
        case 4
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
        case 5
            % PID so just plot emissions and sea ice
            load([pwd filesep 'PID_1_1_1' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'PID_1_1_1' filesep 'state_data' filesep 'all_data_series'])
            
            x1 = find(Param.t==2010);x2 = find(Param.t==2100);x3 = find(Param.t==Param.control_start-1);
            H = figure;
            set(H,'color',[1 1 1])
            plot(Param.t(x1:x2),Data.all_NH_SO2((x1:x2),1,1),'k','linewidth',2)
            xlabel('Year','fontsize',14)
            ylabel('Annual total SO_{2} emissions (\itGt\rm)','fontsize',14)
            title({'SO_2 emissions:'; which_project{project_no}},'fontsize',16)
            fig_file_name = [plot_save_dir filesep [Param.project_name '_emissions']]
            print('-dpng', '-r300', fig_file_name)
            % Sea ice figure
            H2 = figure;
            set(H2,'color',[1 1 1])
            plot(Param.t(x1:x2),Data.the_observations((x1:x2),3)+Param.sea_ice_baseline,'k--','linewidth',2)
            hold on
            % for mc = 1:Param.n_model
            %     plot(Param.t(x1:x2),Data.all_ice_states(x1:x2,1,mc))
            %     %plot(Param.t(time_index+1:end),Data.all_ice_states(time_index,2:end,mc),'r')
            % end
            
            % add the target line
            plot(Param.t(x1:x2),Data.DI_des(x1:x2)+Param.sea_ice_baseline,'k','linewidth',2)
            legend('observed','target','location','northwest')
            xlabel('Year','fontsize',14)
            ylabel({'Minimum sea ice extent';'(\itmillions of square km\rm)'},'fontsize',14)
            title({'Sea ice:'; which_project{project_no}},'fontsize',16)
            fig_file_name = [plot_save_dir filesep [Param.project_name '_sea_ice']]
            print('-dpng', '-r300', fig_file_name)
    end
    % emissions values
    x1 = find(Param.t==2010);x2 = find(Param.t==2100);
    H = figure;
    set(H,'color',[1 1 1])
    P_han = patch([Param.t(x1:x2);flipud(Param.t(x1:x2))],[Data.all_NH_SO2((x1:x2),2,1);flipud(Data.all_NH_SO2((x1:x2),2,10))],[0.5 0.5 0.5])
    set(P_han,'edgecolor','none')
    hold on
    plot(Param.t(x1:x2),Data.all_NH_SO2((x1:x2),1,1),'k','linewidth',2)
    legend('multi-model envelope',char('applied performance','weighted value'),'location','northwest')
    xlabel('Year','fontsize',14)
    ylabel('Annual total SO_{2} emissions (\itGt\rm)','fontsize',14)
    title({'SO_2 emissions:'; which_project{project_no}},'fontsize',16)
    fig_file_name = [plot_save_dir filesep [Param.project_name '_emissions']]
    print('-dpng', '-r300', fig_file_name)
    % Sea ice figure
    H2 = figure;
    set(H2,'color',[1 1 1])
    P_han2 = patch([Param.t(x1:x2);flipud(Param.t(x1:x2))],[Data.all_ice_states((x1:x2),1,1);flipud(Data.all_ice_states((x1:x2),1,10))]+Param.sea_ice_baseline,[0.5 0.5 0.5])
    set(P_han2,'edgecolor',[0.5 0.5 0.5])
    hold on
    plot(Param.t(x1:x2),Data.the_observations((x1:x2),3)+Param.sea_ice_baseline,'k--','linewidth',2)
    hold on
    % for mc = 1:Param.n_model
    %     plot(Param.t(x1:x2),Data.all_ice_states(x1:x2,1,mc))
    %     %plot(Param.t(time_index+1:end),Data.all_ice_states(time_index,2:end,mc),'r')
    % end
    
    % add the target line
    plot(Param.t(x1:x2),Data.DI_des(x1:x2)+Param.sea_ice_baseline,'k','linewidth',2)
    legend(char('multi-model envelope','of the 1 year forecast'),'observed','target','location','northwest')
    xlabel('Year','fontsize',14)
    ylabel({'Minimum sea ice extent';'(\itmillions of square km\rm)'},'fontsize',14)
    title({'Sea ice:'; which_project{project_no}},'fontsize',16)
    fig_file_name = [plot_save_dir filesep [Param.project_name '_sea_ice']]
    print('-dpng', '-r300', fig_file_name)
    
    % temperature figure
    
    H3 = figure;
    set(H3,'color',[1 1 1])
    NH_temp_range = [Data.all_Y_states(1,x1:x2,1) fliplr(Data.all_Y_states(1,x1:x2,10))]+Param.NH_baseline;
    SH_temp_range = [Data.all_Y_states(2,x1:x2,1) fliplr(Data.all_Y_states(2,x1:x2,10))]+Param.SH_baseline;
    NH_adaptive_range_forward = ((Data.all_Y_states(1,x1:x2,1)+Param.offset)'.*Data.all_a_g_NH(x1:x2,1))-Param.offset;
    NH_adaptive_range_backwards = ((Data.all_Y_states(1,x1:x2,10)+Param.offset)'.*Data.all_a_g_NH(x1:x2,10))-Param.offset;
    NH_adaptive_range = [NH_adaptive_range_forward;flipud(NH_adaptive_range_backwards)]+Param.NH_baseline;
    SH_adaptive_range_forward = ((Data.all_Y_states(2,x1:x2,1)+Param.offset)'.*Data.all_a_g_SH(x1:x2,1))-Param.offset;
    SH_adaptive_range_backwards = ((Data.all_Y_states(2,x1:x2,10)+Param.offset)'.*Data.all_a_g_SH(x1:x2,10))-Param.offset;
    SH_adaptive_range = [SH_adaptive_range_forward;flipud(SH_adaptive_range_backwards)]+Param.SH_baseline;
    
    
    P_han3 = patch([Param.t(x1:x2); flipud(Param.t(x1:x2))],SH_temp_range,[0.6 0.6 0.6])
    set(P_han3,'edgecolor',[0.6 0.6 0.6])
    hold on
    P_han4 = patch([Param.t(x1:x2); flipud(Param.t(x1:x2))],NH_temp_range,[0.5 0.5 0.5])
    set(P_han4,'edgecolor',[0.5 0.5 0.5])
    if is_adaptive
        P_han5 = patch([Param.t(x1:x2); flipud(Param.t(x1:x2))],NH_adaptive_range,[0.4 0.4 0.4])
        set(P_han5,'edgecolor',[0.4 0.4 0.4])
        P_han6 = patch([Param.t(x1:x2); flipud(Param.t(x1:x2))],SH_adaptive_range,[0.3 0.3 0.3])
        set(P_han6,'edgecolor',[0.3 0.3 0.3])
    end
    plot(Param.t(x1:x2),Data.the_observations(x1:x2,1)+Param.NH_baseline,'k','linewidth',2)
    plot(Param.t(x1:x2),Data.the_observations(x1:x2,2)+Param.SH_baseline,'k--','linewidth',2)
    if is_adaptive
        legend('NH multi-model envelope','SH multi-model envelope','NH adaptive multi-model envelope',...
            'SH adaptive multi-model envelope','NH observations','SH observations','location','northwest')
    else
        legend('NH multi-model envelope','SH multi-model envelope',...
            'NH observations','SH observations','location','northwest')
    end
    
    xlabel('Year','fontsize',14)
    ylabel({'Temperature';'(^o\itC\rm)'},'fontsize',14)
    title({'Temperature:'; which_project{project_no}},'fontsize',16)
    fig_file_name = [plot_save_dir filesep [Param.project_name '_temperature']]
    print('-dpng', '-r300', fig_file_name)
    
end


% do heteroskedastic comparison
% get an error series
clear targ obs_ice resids resid_var ice_var

for project_no = 1:5
    
    %plot_save_dir = '/Users/davidleedal/Dropbox/Sea_ice_control_papers/policy_as_recursive_process/PRP_figures'
    
    % load each projects data
    switch project_no
        case 1
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 2
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 3
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
        case 4
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
            case 5
            % PID so just plot emissions and sea ice
            load([pwd filesep 'PID_1_1_1' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'PID_1_1_1' filesep 'state_data' filesep 'all_data_series'])
    end
    t_start = 2025;
    targ(:,project_no) = Data.DI_des(find(Param.t == t_start):find(Param.t == 2100));
    obs_ice(:,project_no) = Data.the_observations(find(Param.t == t_start):find(Param.t == 2100),3);
    resids(:,project_no) = targ(:,project_no) - obs_ice(:,project_no);
    resid_var(:,project_no) = heteroskedastic_variance( resids(:,project_no),5 );
    ice_var(:,project_no) = heteroskedastic_variance( obs_ice(:,project_no),5 );
    
end
figure
plot(Param.t(find(Param.t == t_start):find(Param.t == 2100)),resids)
legend('simple multi-model','multi-model with integral','adaptive multi-model with integral','adaptive multi-model','PID')

figure
plot(Param.t(find(Param.t == t_start):find(Param.t == 2100)),resid_var)
legend('simple multi-model','multi-model with integral','adaptive multi-model with integral','adaptive multi-model','PID')
figure
plot(Param.t(find(Param.t == t_start):find(Param.t == 2100)),ice_var)
legend('simple multi-model','multi-model with integral','adaptive multi-model with integral','adaptive multi-model','PID')

% quick plot of weighted SSE
%     lin_style = {'b' 'g' 'r' 'k' 'm' 'c' 'y' 'b:' 'g:' 'r:'};
%     for i = 1:Param.n_model
%     plot(Param.t,squeeze(Data.model_errors(:,5,i)),lin_style{i})
%     hold on
%     end
% legend('1','2','3','4','5','6','7','8','9','10')
% plot(Param.t,Data.all_forcings(:,[1 4 5 6 8 9 10 11]))
% legend('ghg','net TOA','F aero','NH geo','NH v', 'SH v','total NH','total SH')


% comparison plot of various approaches
    H5 = figure;
    set(H5,'color',[1 1 1])
    index_number = [1 3 5 6];
    col_name = {'b' 'g' 'r' 'c'}; 
    for project_no = 1:length(index_number)
    
    %plot_save_dir = '/Users/davidleedal/Dropbox/Sea_ice_control_papers/policy_as_recursive_process/PRP_figures'
    
    % load each projects data
    switch index_number(project_no)
        case 1
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 2
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_no_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 0;
        case 3
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'with_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
        case 4
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'no_IofE_with_parameter_update' filesep 'state_data' filesep 'all_data_series'])
            is_adaptive = 1;
        case 5
            % PID so just plot emissions and sea ice
            load([pwd filesep 'PID_1_1_1' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'PID_1_1_1' filesep 'state_data' filesep 'all_data_series'])
        case 6
            % pure inversion so just plot emissions and sea ice
            load([pwd filesep 'single_model_no_volcanic' filesep 'MagicC_model_parameter_sets'])
            load([pwd filesep 'single_model_no_volcanic' filesep 'state_data' filesep 'all_data_series'])

    end

    
    plot(Param.t(x1:x2),Data.the_observations((x1:x2),3)+Param.sea_ice_baseline,col_name{project_no},'linewidth',2)
    hold on
    % for mc = 1:Param.n_model
    %     plot(Param.t(x1:x2),Data.all_ice_states(x1:x2,1,mc))
    %     %plot(Param.t(time_index+1:end),Data.all_ice_states(time_index,2:end,mc),'r')
    % end
    end
    % overcolour the section of obs pre-control
    plot(Param.t(x1:x3),Data.the_observations((x1:x3),3)+Param.sea_ice_baseline,'m','linewidth',2)

    % add the target line
    plot(Param.t(x1:x2),Data.DI_des(x1:x2)+Param.sea_ice_baseline,'k','linewidth',2)
    %legend(char('multi-model envelope','of the 1 year forecast'),'observed','target','location','northwest')
    xlabel('Year','fontsize',14)
    ylabel({'Minimum sea ice extent';'(\itmillions of square km\rm)'},'fontsize',14)
    %title({'Sea ice:'; which_project{project_no}},'fontsize',16)
    %fig_file_name = [plot_save_dir filesep [Param.project_name '_sea_ice']]
    %print('-dpng', '-r300', fig_file_name)
    
