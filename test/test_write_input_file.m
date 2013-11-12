function test_write_input_file( base_dir, project_name, year, NH_t, SH_t, sea_i, NH_so2, SH_so2,...
    f_ghg, f_sw_down, f_sw_up, net_TOA, NH_v, SH_v)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% input file MUST looks like this:
row_names = {'year' 'NH_temperature' 'SH_temperature' 'minimum_sea_ice_extent' ...
    'NH_SO2_emission' 'SH_SO2_emission' ...
    'f_ghg' 'f_sw_down' 'f_sw_up' 'net_toa' 'NH_volcano_size' 'SH_volcano_size'};
row_values = [year, NH_t, SH_t, sea_i, NH_so2, SH_so2,...
    f_ghg, f_sw_down, f_sw_up, net_TOA, NH_v, SH_v];
% year,2000
% NH_temperature,13
% SH_temperature,12
% minimum_sea_ice_extent,3
% NH_SO2_emission,5
% SH_SO2_emission,0
% f_ghg,3
% f_sw_down,NA
% f_sw_up,NA
% net_toa,NA
% NH_volcano_size,0
% SH_volcano_size,0
fid = fopen([base_dir filesep project_name filesep 'inputs_and_outputs' filesep 'input_file_' num2str(year)], 'a');
for i = 1:12
    fprintf(fid,'%s\n',[row_names{i} ',' num2str(row_values(i))]);
end
fclose(fid)


end

