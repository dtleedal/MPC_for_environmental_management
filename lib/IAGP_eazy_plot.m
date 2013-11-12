function H = IAGP_eazy_plot( base_directory, Data, Param, year_from, year_to, actual_F_aero )
% year_from = 1990;
% year_to = 2059;
start_index = find(Param.t==year_from);
end_index = find(Param.t==year_to);

H = figure;
% NH and SH temperature plot

plot(Param.t(start_index:end_index),Data.the_observations(start_index:end_index,1:2),'linewidth',2)
hold on
% estimate (with data assimilation)
for mc = 1:Param.n_model
    plot(Param.t(start_index:end_index),Data.all_Y_states(1:2,start_index:end_index,mc))
    plot(Param.t(start_index:end_index),((Data.all_Y_states(1,start_index:end_index,mc)+Param.offset)'.*Data.all_a_g_NH(start_index:end_index,mc))-Param.offset,'k')
    plot(Param.t(start_index:end_index),((Data.all_Y_states(2,start_index:end_index,mc)+Param.offset)'.*Data.all_a_g_SH(start_index:end_index,mc))-Param.offset,'k')

end
xlabel('Year')
ylabel('Temperature (perturbation)')
title({'Observed NH temperature (thick blue), model ensemble NH temperature (blue)'; 'Observed SH temperature (thick green), model ensemble SH temperature (green)'})

print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc')

close(H)

% emissions plot
H = figure;
plot(Param.t(start_index:end_index),Data.all_NH_SO2(start_index:end_index,1,1),'linewidth',2)
hold on
% for mc = 1:Param.n_model
%     plot(Param.t(start_index:end),Data.all_NH_SO2(start_index,2:end,mc),'linewidth',1)
% end
xlabel('Year')
ylabel('NH SO2 emissions (annual)')
title('NH SO2 emissions + forecasts')

print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')

close(H)

% ice plot

H = figure;
plot(Param.t(start_index:end_index),Data.the_observations(start_index:end_index,3),'linewidth',2)
hold on
for mc = 1:Param.n_model
    plot(Param.t(start_index:end_index),Data.all_ice_states(start_index:end_index,1,mc))
    %plot(Param.t(time_index+1:end),Data.all_ice_states(time_index,2:end,mc),'r')
end

% add the target line
plot(Param.t(start_index:end_index),Data.DI_des(start_index:end_index),'k','linewidth',2)
xlabel('Year')
ylabel('Min sea ice extent (perturbation)')
title({'Observed sea ice minimum (thick blue), target (thick black),'; 'model ensemble 1-step forecasts (thin blue)'})

print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')


close(H)

H = figure;
plot(Param.t(start_index:end_index),Data.all_forcings(start_index:end_index,5),'b','linewidth',2)
hold on
plot(Param.t(start_index:end_index),actual_F_aero(start_index:end_index),'k','linewidth',2)
xlabel('Year')
ylabel('Forcing (W m^{-2})')
title('Estimated non-ghg anthropogenic forcing (blue), actual (black)')
print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')
close(H)

% adaptive gain NH temp
H = figure;
plot(Param.t(start_index:end_index),Data.all_a_g_NH(start_index:end_index,:))
xlabel('Year')
ylabel('Adaptive gain for NH temperature')
title('NH temperature gain')
print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')
close(H)

% adaptive gain SH temp
H = figure;
plot(Param.t(start_index:end_index),Data.all_a_g_SH(start_index:end_index,:))
xlabel('Year')
ylabel('Adaptive gain for SH temperature')
title('SH temperature gain')
print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')
close(H)
% adaptive gain sea ice
H = figure;
plot(Param.t(start_index:end_index),Data.all_a_g_SI(start_index:end_index,:))
xlabel('Year')
ylabel('Adaptive gain for sea ice')
title('sea ice gain')
print(H,[base_directory filesep Param.project_name filesep 'figures' filesep 'result_figure_eazy'],'-dpsc', '-append')
close(H)


end


