function MC_parameter = IAGP_return_MC_par(n_model,par_bl,par_dist,par_spread)




switch par_dist
    case 1 % uniform percentage spread
        min_val = par_bl*(100-par_spread)/100;
        max_val = par_bl*(100+par_spread)/100;
        MC_parameter = min_val+rand(n_model,1).*(max_val-min_val);
        
    case 2 % normal distribution with std par_spread
        MC_parameter = par_bl + randn(n_model,1)*par_spread;
    otherwise
        error('unknown MC method')
end
MC_parameter(1) = par_bl;

