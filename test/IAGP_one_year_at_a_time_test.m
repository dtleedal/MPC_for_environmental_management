% testing script for one year at a time

clear all
close all
base_directory = [filesep 'Users' filesep 'davidleedal' filesep 'Dropbox' filesep 'IAGP_HadGEM_control' filesep 'AJs_control_model' filesep 'MPC_for_Lawrence'];
base_directory = pwd;
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

IAGP_data_obs = [NHTannual; SHTannual;...
    ((NHTannual+SHTannual)./2); SIminExtended];
plot(IAGP_data_obs')
startYear = 1860; endYear = 2017;
IAGP_test_forcing = RCP45forcings(RCP45forcingsTime>=startYear & RCP45forcingsTime<=2100);
%plot([IAGP_data_obs' IAGP_test_forcing])
emis = 0


% variables for test model
test_model_X = zeros(80,length(t));
test_model_Y = zeros(4,length(t));
test_model_Y_var = zeros(4,length(t));
test_model_P = diag(ones(80,1).*1000);
assimilate = 1;
[test_model_A,test_model_B,test_model_C,test_model_D,test_model_Q,test_model_R] = make_test_model;

for i = startYear:endYear
    array_index = find(t==i);
    forcings = [IAGP_test_forcing(array_index); IAGP_test_forcing(array_index); emis].*2;
    
    % call to main function
    emissions = IAGPperformMPC(i,IAGP_data_obs(:,i-(startYear-1)),...
        emis,IAGP_test_forcing(i-(startYear-1)),base_directory);
    % ---------------------
    
    %disp([num2str(i) ' : ' num2str(mean(emissions))])
    offset_obs = IAGP_data_obs(:,array_index-1) - [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline];
    
    % call to simulation model
    [test_model_X_hat, test_model_Y_hat,test_model_P,test_model_Y_var] = kalmanFilterSim...
        (test_model_X(:,array_index-1),...
        offset_obs,...
        forcings,...
        dt,...
        test_model_A,test_model_B,test_model_C,test_model_D,test_model_Q,...
        test_model_R,test_model_P,...
        assimilate);
    % ------------------------
    test_model_X(:,array_index) = test_model_X_hat;
    test_model_Y(:,array_index) = test_model_Y_hat;
    test_model_Y_var(:,array_index) = diag(test_model_Y_var);

    
end

i = 2018;% edit this for each year
irw_noise = cumsum(randn(253,1).*0.01);
    array_index = find(t==i);
    forcings = [IAGP_test_forcing(array_index); IAGP_test_forcing(array_index); emissions].*2;
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
    off_set_obs = test_model_Y_hat + [NH_baseline;SH_baseline;global_baseline;sea_ice_baseline]+test_model_noise
    disp(['The forcing was: ' num2str(IAGP_test_forcing(array_index))])
    disp(['The emissions where: ' num2str(emissions)])
    IAGP_main_script
    emissions = 0; % change to the value suggested by IAGP_main_script
    
    
%     emissions = IAGPperformMPC(i,off_set_obs,...
%         emissions,IAGP_test_forcing(array_index-1),base_directory);


% testing one year at a time
