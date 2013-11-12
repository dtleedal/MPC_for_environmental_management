% This script reads the RCP4.5 data from excel and calculates the
% relationship between GHG-based forcing and total direct aerosol forcing
%
% the net forcing including volcanic is used to spin up the controller

% the script only reads the subset between 1859 and 2110 as this matches
% the size of the control states

 total_inc_vol = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','B155:B406');
 total_anth = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','E155:E406');
 total_ghg = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','F155:F406');
 total_dir_aer = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','AP155:AP406');
 cloud = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','AW155:AW406');
 year_base = xlsread('RCP45_MIDYEAR_RADFORCING.xlsx','RCP45_MIDYEAR_RADFORCING','A155:A406');
 figure(3)
 plot(year_base,[total_anth total_ghg total_inc_vol total_dir_aer cloud])
 legend('total anth forcing','total GHG forcing','total inc. volcanic','total direct aerosol','cloud albedo effect','location','northwest')
 
 figure(1)
 i2017 = find(year_base==2017);
 sc = scatter(total_ghg(1:i2017),total_anth(1:i2017));
 set(sc,'marker','o','markeredgecolor','k','markerfacecolor','k','sizedata',20)
 xlabel('total GHG forcing')
 ylabel('total anthropogenic forcing')
 title({'Scatter plot of total GHG forcing against total';'anthropogenic forcing for 1859 - 2017';'together with linear trend'})
 
 figure(2)
 sc2 = scatter(total_ghg(1:i2017),total_dir_aer(1:i2017)+cloud(1:i2017))% shows focing due to aerosol is -0.15 of forcing due to ghg
 set(sc2,'marker','o','markeredgecolor','k','markerfacecolor','k','sizedata',20)
  xlabel('total GHG forcing')
 ylabel('total direct aerosol + cloud albedo effect forcing')
 title({'Scatter plot of total GHG forcing against total';' direct aerosol forcing for 1859 - 2017';'together with linear trend'})
set(gca,'ylim',[-1.4 0.1])
%  % summary and results
%  The total ghg forcing is not at zero in 1860 instead it is at ~0.2
%  Therefore we should subtract this baseline from the forcing that
%  Lawrence inputs (which is the value dtermined by converting CO2, CH4 N2O and cfc's into forcings
%  via IPCC TAR 6.2
%  
%  The total ghg forcing is larger than the total anthropogenic forcing mainly due to aerosols
%  A scatter plot and linear fit suggests total anthropogenic = 0.72 * ghg - 0.079
%  Ignoring the small offset this suggests we should:
%  
%  a) subtract a baseline of 0.2 from the ghg-based forcing
%  b) scale this by 0.72 (at least initially) to find the actual forcing to apply to magicC

% simultanious equations from Andy:
% -0.2 = m*0.2 + c
% -1.3 = m*2.9 + c
% m = -0.41
% c = -0.12