function forcing = IAGP_build_forcing_for_respinup(stop_year)
%stop_year = 2039;
load the_forcing_file
GHG = the_forcing(:,1)
load MagicC_model_parameter_sets
Faero_old_parameters = [-0.41 -0.12];
Faero = nan(size(t));
Faero = Faero_old_parameters(1).*GHG + Faero_old_parameters(2);
Faero((find(t==2010)):(find(t==stop_year))) = -0.5;
forcing = GHG + Faero;

