% This script re-spins-up the controller so that it is as if we have begun
% control again at the present year.
%
% In order to spin up we need the following data:
% 1) HadGEM NH temp
% 2) HadGem SH temp
% 3) HadGEM global temp
% 4) HadGEM annual sea ice minimum extent
% 5) The sea ice target
% 6) The GHG emissions
% 7) The Faero forcing
% 8) The SO2 emissions
%
% The complicating issues are:
% 1) The Faero was originally a function of GHG but switches to an observed
% series in 2010
% 2) Should the target series reset to a new transition path?
% 3) Is the control start year hard coded anywhere?
% 
% The approach for tackling this problem is as follows:
% 
% 1) Build a forcing series by hand and use this as the forcing input to IAGPperformMPC
% function. We then need to move the part of IAGPperformMPC 
% that calculates Faero and effective 
% forcing inside the if (t(time_index) >= control_start) block. In this way
% we can make the control models see any forcing we want until control is
% switched on. The process of building the forcing series should be done in a
% function to keep it tidy.
%
% The obvious approach is to add a switch within IAGPperformMPC that
% changes the new Faero forcing calculation to the old one if the year is
% pre 2020.
%
% need the observed NH, SH, Global and Sea ice data
% get all the pertinent parameters
clear all
close all

load MagicC_model_parameter_sets
path_to_data = '/Users/davidleedal/Dropbox/IAGP_HadGEM_control/AJs_control_model/hadgem_diagnostics';
load([path_to_data filesep 'state_data' filesep 'the_observations_file'])
load([path_to_data filesep 'state_data' filesep 'all_SO2_states_file'])
load([path_to_data filesep 'state_data' filesep 'the_forcing_file'])

% we can also get the atmospheric sulphate AOD
% time_base_m = xlsread([path_to_data filesep 'Sulphate_AOD_with_noise_1980_2045.xlsx'],'A3:A794');
% NH_aod_m = xlsread([path_to_data filesep 'Sulphate_AOD_with_noise_1980_2045.xlsx'],'B3:B794');
% SH_aod_m = xlsread([path_to_data filesep 'Sulphate_AOD_with_noise_1980_2045.xlsx'],'C3:C794');
% time_base = reshape(time_base_m,12,size(NH_aod_m,1)./12);time_base = round(time_base(1,:));
% NH_aod = mean(reshape(NH_aod_m,12,size(NH_aod_m,1)./12));
% SH_aod = mean(reshape(SH_aod_m,12,size(NH_aod_m,1)./12));
% plot(time_base,NH_aod)
% NH_aod_sub = NH_aod(find(time_base==2010):find(time_base==2044))

% we also need a volcano series
% I'm going to build it here
volcano = zeros(size(t));
% Pinatubo
volcano(find(t==1980))=5;
volcano(find(t==1982))=5;
volcano(find(t==1991))=18;
volcano(find(t==2038))=6;
volcano(find(t==2044))=9;






























% testing script for IAGP walk through
%base_directory = [filesep 'Users' filesep 'davidleedal' filesep 'Dropbox' filesep 'IAGP_HadGEM_control' filesep 'AJs_control_model' filesep 'MPC_for_Lawrence'];
base_directory = pwd
IAGP_model_setup(base_directory)
IAGP_data_setup_ONE_OFF(base_directory)


NH_annual = the_observations(1,:) + NH_baseline
SH_annual = the_observations(2,:) + SH_baseline;
global_annual = the_observations(3,:) + mean([NH_baseline SH_baseline]);
sea_ice_min = the_observations(4,:) + sea_ice_baseline

% The state data files all start with zero initial conditions,
% which coincides with the year 1859. Therefore the timeseries
% that start in 1860 need a zero appending at the begining
IAGP_data_obs = [[0;0;0;0] [NH_annual; SH_annual;...
    global_annual; sea_ice_min]];
%plot(IAGP_data_obs')
startYear = 1860; endYear = 2045;
emis = all_SO2(:,1);

for i = startYear:endYear
    array_index = find(t==i);
    
    
    % call to main function
    emissions = IAGPperformMPC(i,IAGP_data_obs(:,array_index),...
        emis(array_index),the_forcing(array_index,1),volcano(array_index),base_directory);
    % ---------------------
    
    
    
    
    
end

return
% irw_noise = zeros(253,1);
% irw_noise = [zeros(159,1); cumsum(randn(253-159,1).*0.002)];
% 
% % testing one year at a time
% 
% 
% 
% emissions = 5;
% for i = endYear+1:2028%2099
%     for i = 2068:2070
%     array_index = find(t==i);
%     forcings = [reduced_IAGP_test_forcing(array_index); reduced_IAGP_test_forcing(array_index); emissions].*2;
%     [test_model_X_hat, test_model_Y_hat,test_model_P,test_model_Y_var] = kalmanFilterSim...
%         (test_model_X(:,array_index-1),...
%         test_model_Y(:,array_index-1),...
%         forcings+irw_noise(array_index),...
%         dt,...
%         test_model_A,test_model_B,test_model_C,test_model_D,...
%         test_model_Q, test_model_R,test_model_P,...
%         0);
%     test_model_X(:,array_index) = test_model_X_hat;
%     test_model_Y(:,array_index) = test_model_Y_hat;
%     test_model_Y_var(:,array_index) = diag(test_model_Y_var);
%     test_model_noise = (2.*test_model_R)*randn(4,1);
%     off_set_obs = test_model_Y_hat + [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline]+test_model_noise;
%     emissions = IAGPperformMPC(i,off_set_obs,...
%         emissions,IAGP_test_forcing(array_index),base_directory);
% end
% 
% % load the stored states
% 
% 
load([base_directory filesep 'state_data' filesep 'all_Y_states_file'])
load([base_directory filesep 'state_data' filesep 'all_SO2_states_file'])
load([base_directory filesep 'state_data' filesep 'DI_des_file'])%           DI_des
load([base_directory filesep 'state_data' filesep 'all_Y_states_SE_file'])
load([base_directory filesep 'state_data' filesep 'the_forcing_file'])
load([base_directory filesep 'state_data' filesep 'the_observations_file'])
%bar(reshape(all_Y_states_SE(1,find(t==2017),:),n_model,1))
%bar(reshape(all_Y_states_SE(4,find(t==2017),:),n_model,1))
% 
% % t is the time index up to 2100 + N
% % define the start and end points for plotting
 t1860 = find(t==1860); t1990 = find(t==1990); t2017 = find(t==2017);t2044 = find(t==2044);t2099 = find(t==2099);
% 
% 
% 
which_n = 1;
which_state = 4;
% 
 plot(t(t1860:t2017)',IAGP_data_obs(which_state,t1860:t2017)')
hold on
plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)')
plot(t(t1860:t2099)',test_model_Y(which_state,t1860:t2099)','k','linewidth',2)
% 
plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)'+2.*sqrt(all_Y_states_SE(which_state,t1860:t2099,which_n)'),'b:')
plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)'-2.*sqrt(all_Y_states_SE(which_state,t1860:t2099,which_n)'),'b:')
plot(t(t1860:t2099)',DI_des(t1860:t2099)','k')
plot(t(t1860:t2099)',the_observations(which_state,t1860:t2099)')

plot(t(t1860:t2099)',all_SO2(t1860:t2099,1,which_n)','r')

scatter(the_observations(1,t1990:t2044),the_observations(4,t1990:t2044))
hold on
scatter(all_Y_states(1,t1990:t2044,which_n),all_Y_states(4,t1990:t2044,which_n),'ro')

% %nextForcing = GasConversionFunction(CO2_ppm,CH4_ppm,N2O_ppm,CFC_11_ppb,CFC_12_ppb)
% % convert gas to forcing here?
% nextForcing = GasConversionFunction(350,1.8,0.4,0,0)
% 
% % should we subtract baselines here?
% NHtemp = 14.38 - baselineNHT;
% SHtemp = 14.3 - baselineSHT;
% Gtemp = ((NHtemp + SHtemp)/2) - baselineGT;
% SI = 3.1 - baselineSImin;
% theYear = 2020;
% emissions = IAGPperformMPC(theYear,[NHtemp; SHtemp; Gtemp; SI],26,nextForcing,this_directory);
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% centres = 0:0.25:20;
% emissions_frequencies = zeros(length(centres),t2100-(t2020-1));
% for i = t2020:t2100
%     emissions_frequencies(:,i-(t2020-1)) = hist(reshape(all_SO2(i,2,:),1,n_model),centres);
% end
% emissions_frequencies(emissions_frequencies>50) = 50;
% fhan = figure
% set(fhan,'color',[1 1 1])
% emisHan = plot3(t2020:t2100,all_SO2(t2020:t2100,2,1),ones(size(all_SO2(t2020:t2100,2,1))).*60)
% set(emisHan,'color','k','linewidth',2)
% hold on
% surfHan = surface(t2020:t2100,centres,emissions_frequencies);
% set(surfHan,'linestyle','none','facecolor','interp')
% set(gca,'xlim',[2020 2100],'ylim',[0 20],'layer','top')
% load colormapForHeatplot
% colormap(colormapForHeatplot)
% cbHan = colorbar('east')
% set(cbHan,'position',[0.8 0.15 0.08 0.3],'Ytick',[0 10 20 30 40 50],'YtickLabel'...
%     ,[0 10 20 30 40 50]./n_model)
% grid on
% 
% xlabel('year')
% ylabel('SO_2 emissions (\it{GtC/year}\rm)')
% title('Probabilistic SO_2 emissions to achieve sea ice stabilization (m=500)')
% legHan = legend('mean','location','northwest')
% 
% 
% 
% 
% 
% 
% % % sea ice extent
% % centres = -0.22:0.003:0.06;
% % sea_ice_frequencies = zeros(length(centres),190-121);
% % for i = 122:190
% %     sea_ice_frequencies(:,i-121) = hist(reshape(all_Y_states(end,i,:),1,500),centres);
% % end
% % %sea_ice_frequencies(sea_ice_frequencies>50) = 50;
% % fhan2 = figure
% % set(fhan2,'color',[1 1 1])
% % targetHan = plot3(t(122):t(190),DI_des(122:190),ones(size(DI_des(122:190))).*100)
% % set(targetHan,'color','k','linewidth',2)
% % surfHan2 = surface(t(122):t(190),centres,sea_ice_frequencies);
% % set(surfHan2,'linestyle','none','facecolor','interp')
% % set(gca,'xlim',[2020 2090],'ylim',[-0.22 0.06],'layer','top')
% % load colormapForHeatplot
% % colormap(colormapForHeatplot)
% % cbHan2 = colorbar('east')
% % set(cbHan2,'position',[0.8 0.15 0.08 0.3],'Ytick',[0 20 40 60 80],'YtickLabel'...
% %     ,[0 20 40 60 80]./500)
% % grid on
% %
% % xlabel('year')
% % ylabel('\Deltasea ice extent (\itMkm^2\rm)')
% % title('Probabilistic sea ice stabilization (m=500)')
% % legHan = legend('target','location','northwest')
% %
