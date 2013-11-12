function [year,this_obs_original,NH_SO2_emission,SH_SO2_emission,f_ghg,f_sw_down,f_sw_up,net_toa,NH_volcano_size,SH_volcano_size]...
    = IAGP_read_input_file(year_from_cmd,base_directory,project_name)



% input file MUST looks like this:

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

disp(['Loading data for year: ' num2str(year_from_cmd)])
i = 1;
file_to_open = [base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'input_file_' num2str(year_from_cmd)];
try
    fid=fopen(file_to_open,'r');
    while ~feof(fid)
        tline = fgetl(fid);
        try
            C = regexp(tline,',','split');
            if ~isempty(find(strcmp(C{2},{'NAN', 'NaN', 'nan', 'NA', 'na', 'Na', '-9999', '', ' ', '-'}))==1)
                disp(['Warning: detected a missing value for ' C{1} '. Registering this as NaN'])
                C{2} = 'NaN';
            end
            eval([C{1} ' = ' C{2} ';']);
        catch
            error(['Could not convert input_file line number ' num2str(i) ' into matlab variable'])
        end
        i = i + 1;
    end
    fclose(fid);
catch
    error(sprintf(['Could not open: ' file_to_open '\n'...
        'possible reasons:\n(1) year from command line greater than year of input_file\n'...
        '(2) problem with the current working directory --you shoud be in <path to>/MPC_for_Ben']));
end
% put the temperatures and sea ice into a vector
this_obs_original = [NH_temperature;SH_temperature;minimum_sea_ice_extent];


