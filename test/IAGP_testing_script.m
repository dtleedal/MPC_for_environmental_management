% testing script for IAGP walk through
clear all
close all
%base_directory = [filesep 'Users' filesep 'davidleedal' filesep 'Dropbox' filesep 'IAGP_HadGEM_control' filesep 'AJs_control_model' filesep 'MPC_for_Lawrence'];
base_directory = pwd
IAGP_model_setup(base_directory)
IAGP_data_setup_ONE_OFF(base_directory)
load RCP45forcings
load Hadley_data_for_andy_v1_3
load([base_directory filesep 'MagicC_model_parameter_sets'])
% the RCP8.5 data begins in 1765
% the HadGEM2 data begins 1860
% so the first RCP8.5 value matching the HadGEM2 data is index 96
% and 2017 (where HadGEM2 data ends) is index 253
% put the HadGEM2 data into offset form and place into an observations
% matrix
%plot([NHTannual' SHTannual' (NHTannual'+SHTannual')./2])
% baselineNHT = mean(NHTannual(1:40));
% baselineSHT = mean(SHTannual(1:40));
% baselineGT = mean((NHTannual(1:40)+SHTannual(1:40))./2);

% baselineSImin = 5.5;

SIminExtended = (-2.*(NHTannual-NH_baseline))+sea_ice_baseline;
SIminExtended(end-(length(SImin)-1):end) = (SImin);
%plot(SIminExtended)

% The state data files all start with zero initial conditions,
% which coincides with the year 1859. Therefore the timeseries
% that start in 1860 need a zero appending at the begining
IAGP_data_obs = [[0;0;0;0] [NHTannual; SHTannual;...
    ((NHTannual+SHTannual)./2); SIminExtended]];
plot(IAGP_data_obs')
startYear = 1860; endYear = 2017;
IAGP_test_forcing = RCP45forcings(RCP45forcingsTime>=(startYear-1) & RCP45forcingsTime<=2110);
reduced_IAGP_test_forcing = IAGP_test_forcing + forcing_conversion(1).*IAGP_test_forcing + forcing_conversion(2);
%IAGP_test_forcing = IAGP_test_forcing - mean(IAGP_test_forcing(1:40));
%plot([IAGP_data_obs' IAGP_test_forcing])
emis = 0
[test_model_A,test_model_B,test_model_C,test_model_D,test_model_Q,test_model_R] = make_test_model;

% variables for test model
test_model_X = zeros(80,length(t));
test_model_Y = zeros(4,length(t));
test_model_Y_var = zeros(4,length(t));
test_model_P = diag(ones(80,1).*1000);
assimilate = 0;

% volcanoes
test_volcanoes = zeros(length(t),1);
test_volcanoes(find(t==2038)) = 6;
test_volcanoes(find(t==2044)) = 9;


% then overwrite test model with one of the 10 control models
% test_model_A = A(:,:,1);
% test_model_B = B(:,:,1);
% test_model_C = C(:,:,1);
% test_model_D = D(:,:,1);

% then make synthetic HadGEM observations using the test model as a stand
% in
% irw_noise = cumsum(randn(159,1).*0.01);
% for i = startYear:endYear
%     array_index = find(t==i);
%     forcings = [reduced_IAGP_test_forcing(array_index); reduced_IAGP_test_forcing(array_index); 0].*2;
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
%     test_model_noise = (1.5.*test_model_R)*randn(4,1);
%     IAGP_data_obs(:,array_index) = test_model_Y_hat + [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline]+test_model_noise;
% end
% 

for i = startYear:endYear
    array_index = find(t==i);
    
    forcings = [reduced_IAGP_test_forcing(array_index); reduced_IAGP_test_forcing(array_index); emis].*2;
    offset_obs = IAGP_data_obs(:,array_index) - [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline];
    [test_model_X_hat, test_model_Y_hat,test_model_P,test_model_Y_var] = kalmanFilterSim...
        (test_model_X(:,array_index-1),...
        offset_obs,...
        forcings,...
        dt,...
        test_model_A,test_model_B,test_model_C,test_model_D,test_model_Q,...
        test_model_R,test_model_P,...
        assimilate);
    
    volcano_size = test_volcanoes(array_index);
    % call to main function
    emissions = IAGPperformMPC(i,IAGP_data_obs(:,array_index),...
        emis,IAGP_test_forcing(array_index),volcano_size,base_directory);
    % ---------------------
    
    
    
    % call to simulation model
    % ------------------------
    test_model_X(:,array_index) = test_model_X_hat;
    test_model_Y(:,array_index) = test_model_Y_hat;
    test_model_Y_var(:,array_index) = diag(test_model_Y_var);
    
    
end
irw_noise = zeros(253,1);
irw_noise = [zeros(159,1); cumsum(randn(253-159,1).*0.002)];

% testing one year at a time



%emissions = 5;
for i = endYear+1:2045%2099
    for i = 2037:2045
    array_index = find(t==i);
    forcings = [reduced_IAGP_test_forcing(array_index); reduced_IAGP_test_forcing(array_index); emissions].*2;
    [test_model_X_hat, test_model_Y_hat,test_model_P,test_model_Y_var] = kalmanFilterSim...
        (test_model_X(:,array_index-1),...
        test_model_Y(:,array_index-1),...
        forcings+irw_noise(array_index),...
        dt,...
        test_model_A,test_model_B,test_model_C,test_model_D,...
        test_model_Q, test_model_R,test_model_P,...
        0);
    test_model_X(:,array_index) = test_model_X_hat;
    test_model_Y(:,array_index) = test_model_Y_hat;
    test_model_Y_var(:,array_index) = diag(test_model_Y_var);
    test_model_noise = (2.*test_model_R)*randn(4,1);
    off_set_obs = test_model_Y_hat + [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline]+test_model_noise;
    volcano_size = test_volcanoes(array_index);
    emissions = IAGPperformMPC(i,off_set_obs,...
        emissions,IAGP_test_forcing(array_index),volcano_size,base_directory);
end

% load the stored states


load([base_directory filesep 'state_data' filesep 'all_Y_states_file'])
load([base_directory filesep 'state_data' filesep 'all_SO2_states_file'])
load([base_directory filesep 'state_data' filesep 'DI_des_file'])%           DI_des
load([base_directory filesep 'state_data' filesep 'all_Y_states_SE_file'])
load([base_directory filesep 'state_data' filesep 'the_forcing_file'])
%bar(reshape(all_Y_states_SE(1,find(t==2017),:),n_model,1))
%bar(reshape(all_Y_states_SE(4,find(t==2017),:),n_model,1))

% t is the time index up to 2100 + N
% define the start and end points for plotting
t1860 = find(t==1860); t2017 = find(t==2017);t2099 = find(t==2099);



which_n = 1;
which_state = 1;

plot(t(t1860:t2017)',IAGP_data_obs(which_state,t1860:t2017)')
hold on
plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)')
plot(t(t1860:t2099)',test_model_Y(which_state,t1860:t2099)','k','linewidth',2)

plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)'+2.*sqrt(all_Y_states_SE(which_state,t1860:t2099,which_n)'),'b:')
plot(t(t1860:t2099)',all_Y_states(which_state,t1860:t2099,which_n)'-2.*sqrt(all_Y_states_SE(which_state,t1860:t2099,which_n)'),'b:')
plot(t(t1860:t2099)',DI_des(t1860:t2099)','k')
plot(t(t1860:t2099)',all_SO2(t1860:t2099,1,which_n)','r')
%nextForcing = GasConversionFunction(CO2_ppm,CH4_ppm,N2O_ppm,CFC_11_ppb,CFC_12_ppb)
% convert gas to forcing here?
nextForcing = GasConversionFunction(350,1.8,0.4,0,0)

% should we subtract baselines here?
NHtemp = 14.38 - baselineNHT;
SHtemp = 14.3 - baselineSHT;
Gtemp = ((NHtemp + SHtemp)/2) - baselineGT;
SI = 3.1 - baselineSImin;
theYear = 2020;
emissions = IAGPperformMPC(theYear,[NHtemp; SHtemp; Gtemp; SI],26,nextForcing,this_directory);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centres = 0:0.25:20;
emissions_frequencies = zeros(length(centres),t2100-(t2020-1));
for i = t2020:t2100
    emissions_frequencies(:,i-(t2020-1)) = hist(reshape(all_SO2(i,2,:),1,n_model),centres);
end
emissions_frequencies(emissions_frequencies>50) = 50;
fhan = figure
set(fhan,'color',[1 1 1])
emisHan = plot3(t2020:t2100,all_SO2(t2020:t2100,2,1),ones(size(all_SO2(t2020:t2100,2,1))).*60)
set(emisHan,'color','k','linewidth',2)
hold on
surfHan = surface(t2020:t2100,centres,emissions_frequencies);
set(surfHan,'linestyle','none','facecolor','interp')
set(gca,'xlim',[2020 2100],'ylim',[0 20],'layer','top')
load colormapForHeatplot
colormap(colormapForHeatplot)
cbHan = colorbar('east')
set(cbHan,'position',[0.8 0.15 0.08 0.3],'Ytick',[0 10 20 30 40 50],'YtickLabel'...
    ,[0 10 20 30 40 50]./n_model)
grid on

xlabel('year')
ylabel('SO_2 emissions (\it{GtC/year}\rm)')
title('Probabilistic SO_2 emissions to achieve sea ice stabilization (m=500)')
legHan = legend('mean','location','northwest')






% % sea ice extent
% centres = -0.22:0.003:0.06;
% sea_ice_frequencies = zeros(length(centres),190-121);
% for i = 122:190
%     sea_ice_frequencies(:,i-121) = hist(reshape(all_Y_states(end,i,:),1,500),centres);
% end
% %sea_ice_frequencies(sea_ice_frequencies>50) = 50;
% fhan2 = figure
% set(fhan2,'color',[1 1 1])
% targetHan = plot3(t(122):t(190),DI_des(122:190),ones(size(DI_des(122:190))).*100)
% set(targetHan,'color','k','linewidth',2)
% surfHan2 = surface(t(122):t(190),centres,sea_ice_frequencies);
% set(surfHan2,'linestyle','none','facecolor','interp')
% set(gca,'xlim',[2020 2090],'ylim',[-0.22 0.06],'layer','top')
% load colormapForHeatplot
% colormap(colormapForHeatplot)
% cbHan2 = colorbar('east')
% set(cbHan2,'position',[0.8 0.15 0.08 0.3],'Ytick',[0 20 40 60 80],'YtickLabel'...
%     ,[0 20 40 60 80]./500)
% grid on
%
% xlabel('year')
% ylabel('\Deltasea ice extent (\itMkm^2\rm)')
% title('Probabilistic sea ice stabilization (m=500)')
% legHan = legend('target','location','northwest')
%
