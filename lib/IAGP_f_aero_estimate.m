function [H, f_aero] = IAGP_f_aero_estimate( year,...
    all_forcings,...
    the_observations,...
    all_SO2,...
    two_x_CO2_eq_T,...
    a_g_NH,...
    a_g_SH,...
    offset,...
    SO_initial,...
    phi_bl,...
    do_plot_flag)
% a digitized version of the SO2 to forcing coefficient prob. dens. func. from Andy's CMIP
% analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
%                     README                     %
%                                                %
%   AJ: just load Faero_workspace and you will   %
%   have all the variables generated by this     %
%   script                                       %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all close all
% end_year = 2069; % the last year
% end_row = '62'; % the excel final row
% file_name = 'toa_data_with_noise.xlsx';
% sheet_name = '2010-2069 (8 July)';
% digitizedData = [
%     -0.4            0
%     -0.4112         0
%     -0.3824    0.0106
%     -0.3616    0.0459
%     -0.3392    0.1378
%     -0.3248    0.2367
%     -0.3120    0.3216
%     -0.3040    0.4099
%     -0.2896    0.5053
%     -0.2832    0.5548
%     -0.2640    0.5866
%     -0.2480    0.6714
%     -0.2384    0.6961
%     -0.2272    0.8339
%     -0.2192    0.8940
%     -0.2080    0.9399
%     -0.2032    0.9859
%     -0.1920    0.9399
%     -0.1840    0.8763
%     -0.1744    0.7703
%     -0.1648    0.5654
%     -0.1584    0.5406
%     -0.1536    0.4276
%     -0.1472    0.3322
%     -0.1360    0.2120
%     -0.1328    0.1661
%     -0.1248    0.1131
%     -0.1136    0.0565
%     -0.0976    0.0177
%     -0.0800    0.0071
%     -0.0512    0.0035
%     -0.0400         0
%     -0.03           0
%     ];
% digitizedData = digitizedData + 0.04;% change due to HadGEM2 value
% % parameterize the empirical curve with a spline
% PP = pchip(digitizedData(:,1),digitizedData(:,2));
% % make a CDF
% xvals = -0.4:0.005:-0.04;
% yvals = ppval(PP,xvals);
% yvals = yvals./sum(yvals);
% plot(xvals,yvals);
% title('non-parametric pdf for SO2 coefficient')

% % parameterize CDF with spline
% [C,C1]=unique(cumsum(yvals));
% PP2 = pchip(C,xvals(C1));
% figure
% plot(0:0.01:1,ppval(PP2,0:0.01:1),'r')
% title('inverse CDF for easy sampling','fontsize',16)



% The old net TOA forcing

% N = [1.36; -1.28; 1.91; -2.49; -0.05; -0.46; 1.83; -0.14; 2.43; 0.21; 1.16; 2.00; -1.23; ...
%     1.05; 3.42; 3.06; 1.71; 1.92; 1.96; 5.00; 6.02; 2.83; 1.69; 4.54; 1.00; -5.24; -1.77];

% the new net TOA forcing
% N = [0.84; 1.19; 1.33; 1.20; 0.81; 0.88; 1.08; 1.08; 1.09; 0.77; 1.01;...
% 0.99; 0.67; 1.12; 1.16; 1.84; 1.52; 1.05; 0.99; 1.50; 0.66; 1.30; 0.67;...
% 1.16; 1.30; 1.03; 0.93];
Fghg = all_forcings(:,1); % GHG forcing
Fv_NH = all_forcings(:,8); % NH volcanic forcing
Fv_SH = all_forcings(:,9); % SH volcanic forcing
Fv = Fv_NH + Fv_SH;
SW_in = all_forcings(:,2);   % SW down (into Earth)
SW_out = all_forcings(:,3);   % SW up (out to space)
N = all_forcings(:,4);    % net TOA forcing

% albedo = SW_out./SW_in;

% Fghg = [2.569; 2.602; 2.635; 2.668; 2.702; 2.735; 2.768; 2.802; 2.835; 2.868; ...
%     2.901; 2.934; 2.966; 2.999; 3.032; 3.065; 3.098; 3.132; 3.165; 3.199; 3.232; ...
%     3.265; 3.298; 3.330; 3.362; 3.393; 3.425];


% T = [14.1690; 14.0730; 14.2610; 14.2080; 14.2320; 14.1750; 14.3300; ...
%     14.3820; 14.3750; 14.2600; 14.3220; 14.4460; 14.2640; 14.3740;...
%     14.3900; 14.4290; 14.5260; 14.6510; 14.6340; 14.6060;...
%     14.6160; 14.5780; 14.6830; 14.6300; 14.6140; 14.6550; 14.5830];

%T = xlsread(file_name, sheet_name, ['G3:G' end_row '']);
%DT = mean(the_observations(:,1:2),2);
% [m_dim,n_dim,p_dim] = size(model_temperatures);
% mean_model_temperatures = zeros(n_dim,p_dim);
% for j = 1:p_dim
%     mean_model_temperatures(:,j) = mean(model_temperatures(:,:,j);
% end
% baseline global temperature


% Emis = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0.0374;...
%     0.2287; 0.5186; 0.9083; 1.2144; 1.6413; 1.9672; 2.2473;...
%     2.5946; 3.0651; 3.2788];

Emis = all_SO2;
% the time base 2010 to 2036
%year = xlsread(file_name, sheet_name, ['A3:A' end_row '']);

% 10000 draws from SO2 coefficient distribution
%num_of_MC = 10;%000;
%MCFSO2 = ppval(PP2,rand(num_of_MC,1));
% figure
% hist(MCFSO2,40)

% define the 5 values of Lambda used by the 10 MAGICC models
% the numbers from 4.3193 to 2.0436 are the 2xCO2 eq. sensitivity
L = 3.74./two_x_CO2_eq_T;

% pre-assign memory for MC ensemble
% Faero_ensemble = zeros(num_of_MC,length(DT));
% Faero_medians = zeros(length(two_x_CO2_eq_T),length(DT));
Faero_ensemble = zeros(length(L),length(Fghg));
% TODO for testing set MCFSO2 to -0.16
%MCFSO2 = -0.16;

% convert the SSE performance of each model into an exponetially weighted
% measure to emphasise the contribution from well-fitting models
% [aa,bb] = sort(SSE_performance);
% SSE_performance(bb(end)) = 1;
% SSE_performance(bb(1:end-1)) = 0;

% SSE_performance = SSE_performance./max(SSE_performance);
% exponential_scores = SSE_performance.^3;
% sum_of_exponential_scores = sum(exponential_scores);
% exponential_scores = exponential_scores./sum_of_exponential_scores;

% Then make one Monte Carlo ensemble per MAGICC model of estimates for Faero
%net_TOA = dir_aero(time_ind) + f_ghg + phi_bl_test.*NH_so2 - lambda_delta_T - delta_SI(time_ind);
for i = 1:length(two_x_CO2_eq_T)
    DT_NH = ((the_observations(:,1)+offset).*(1./a_g_NH(:,i)))-10;
    DT_SH = ((the_observations(:,2)+offset).*(1./a_g_SH(:,i)))-10;
    DT = mean([DT_NH DT_SH],2);
    for j = 1:length(DT)
        Faero_ensemble(i,j) = (N(j) - Fghg(j) - Fv(j) -(Emis(j).*phi_bl)./2 + L(i).*DT(j) + (SW_out(j) - SO_initial));%(albedo(j)-.29)*SI(j);
    end
    Faero = sum(Faero_ensemble,1)./length(two_x_CO2_eq_T);
end

% perform IRWSM to fit line through median values.
Faero_smooth = irwsm(Faero,1,1e-2);
f_aero = Faero_smooth(end);
if f_aero <-1.5   
    disp(['WARNING: aerosol forcing is: ' num2str(f_aero) '. This value is too low setting to -1.5'])
    f_aero = -1.5;
elseif f_aero > 0
        disp(['WARNING: aerosol forcing is: ' num2str(f_aero) '. This value is too high setting to 0'])
        f_aero = 0;
end 


% if time_index == 202
%     f_aero = 0;
% end


if do_plot_flag
    H = figure;
    plot(year,Faero,'k-',year,Faero_smooth,'r-')
    title(['Ensemble median for aerosol forcing (' num2str(year(end)) ')'],'fontsize',14)
    xlabel('Year','fontsize',12)
    ylabel('Forcing (\itW m^{-2}\rm)','fontsize',12)
    legend('data','IRWSM')
else
    H = 0; % need a dummy variable for H
end % end of do plot

% % for testing
% f_aero = -0.5;
% H = 0;

% plot the results
% figure
% subplot(2,1,1)
% boxplot(Faero2,'label',year)
% xlabel('year','fontsize',16)
% ylabel('Faero W/m^2','fontsize',16)
% title('Box and whisker plot for estimate of Faero where   \lambda = 3.7/3.3446','fontsize',16)
% subplot(2,1,2)
% plot(year,median(Faero1),'ro',year,median(Faero2),'go',year,median(Faero3),'bo',...
%     year,median(Faero4),'ko',year,median(Faero5),'mo')
% hold on
% line_handle = line([2009.5 end_year+0.5], [-0.5 -0.5]);
% set(line_handle,'linestyle','--','color','k')
% axis([2009.5 end_year+0.5 -8 8])
% legend('\lambda = 3.7/4.32','\lambda = 3.7/3.34','\lambda = 3.7/2.75',...
%     '\lambda = 3.7/2.34','\lambda = 3.7/2.04')
% xlabel('year','fontsize',16)
% ylabel('Faero W/m^2','fontsize',16)
% title('Estimate of median Faero for all 5 climate sensitivity values','fontsize',16)

% then, at this stage, I'm choosing an approximate constant for this fairly stationary series


% example for one version of MAGICC for one simulation year (2036)
% figure
% hist(Faero3(:,end),40)
% title('Example distribution of Faero for 2036 using Model 3')

