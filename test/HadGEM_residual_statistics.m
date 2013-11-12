% estimating noise covariance parameters for KF
load Hadley_data_for_andy_v1_3
% find the data trend
NHTannnualS = irwsm(NHTannual,1,1e-2);
SHTannnualS = irwsm(SHTannual,1,1e-2);
SIminS = irwsm(SImin,1,1e-2);
GTannual = (NHTannual + SHTannual)./2;
GTannualS = irwsm(GTannual,1,1e-2);

% check the trend looks right
plot([NHLannual' NHLannnualS])

residN = NHTannual' - NHTannnualS;
residS = SHTannual' - SHTannnualS;
residI = SImin' - SIminS;
residGT = GTannual' - GTannualS;

qqplot([residN residS residGT])
legend('NH','SH','G')
qqplot([residI])


disp(['Mean residual for NH: ' num2str(mean(residN)) ', variance: ' num2str(var(residN))])
disp(['Mean residual for SH: ' num2str(mean(residS)) ', variance: ' num2str(var(residS))])
disp(['Mean residual for global: ' num2str(mean(residGT)) ', variance: ' num2str(var(residGT))])
disp(['Mean residual for min ice: ' num2str(mean(residI)) ', variance: ' num2str(var(residI))])

% Mean residual for NH: -9.2212e-08, variance: 0.01713
% Mean residual for SH: -4.0088e-06, variance: 0.010838
% Mean residual for global: -2.0505e-06, variance: 0.010025
% Mean residual for min ice: 0.00012007, variance: 0.14487

hist([residN residS residGT],20)
legend('NH','SH','G')

% KS test to check if normal % they all pass at 5% significance
kstest((residN-mean(residN))/std(residN))
kstest((residS-mean(residS))/std(residS))
kstest((residGT-mean(residGT))/std(residGT))
kstest((residI-mean(residI))/std(residI))