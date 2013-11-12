function err = IAGP_get_ice_err(g,i,sea_ice_gain_nodes,which_node,year_series,all_temp,sea_ice)

    
    PP_gain_hat = pchip([sea_ice_gain_nodes(1:which_node-1,1); i],[sea_ice_gain_nodes(1:which_node-1,2); g]);

err = sea_ice - all_temp.*ppval(PP_gain_hat,year_series);
end
