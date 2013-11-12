function fitted_DI_des = IAGP_build_fitted_DI_des(all_ice_states,...
        t,...
        control_start,...
        ice_stabilisation_level,...
        ice_stabilisation_year,...
        sea_ice_baseline)
smooth_end = find(t==control_start);
stabilisation_year_index = find(t==ice_stabilisation_year);
end_ice_level = ice_stabilisation_level - sea_ice_baseline;
[sm_ice, Dsm_ice] = irwsm(all_ice_states(1:smooth_end), 1, 1e-2);  
Y = [Dsm_ice(end) sm_ice(end) mean([sm_ice(end) end_ice_level]) end_ice_level 0];
X = [1 round((ice_stabilisation_year-control_start)/2) ice_stabilisation_year-control_start];
XX = 1:ice_stabilisation_year-control_start;
YY = spline(X,Y,XX);
fitted_DI_des = nan(size(t));
fitted_DI_des(smooth_end:stabilisation_year_index-1,1) = YY';
fitted_DI_des(stabilisation_year_index:end,1) = end_ice_level;

end

