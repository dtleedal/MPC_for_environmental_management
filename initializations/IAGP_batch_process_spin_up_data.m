function [Data, Param] = IAGP_batch_process_spin_up_data(base_directory, project_name, Data, Param)


% IAGPperformMPC(year_from_cmd, base_directory, project_name, Data, Param, spin_up_year, save_data_flag, write_SO2_flag)
%% Load the spin up data to get the start date and check spin up exists
try
    load([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'spin_up_data.mat'])
    
catch
    err_msg = 'File load failed. Can''t load the spin up data file. There should\nbe a .mat file called ''spin_up_data.mat'' containing a matrix called ''spin_up_data''\nsaved in the project''s ''inputs_and_outputs'' directory. Please check.';
    error('errorMsg:converted',err_msg)
end

%% make the first input_file

%% main processing loop

row_index = 1;
for the_year = Param.start_year:Param.control_start-1
    spin_up_year.year = spin_up_data(row_index,1);
    spin_up_year.this_obs_original = [spin_up_data(row_index,2) spin_up_data(row_index,3) spin_up_data(row_index,4)];
    spin_up_year.NH_SO2_emission = spin_up_data(row_index,5);
    spin_up_year.SH_SO2_emission = spin_up_data(row_index,6);
    spin_up_year.f_ghg = spin_up_data(row_index,7);
    spin_up_year.f_sw_down = spin_up_data(row_index,8);
    spin_up_year.f_sw_up = spin_up_data(row_index,9);
    spin_up_year.net_toa = spin_up_data(row_index,10);
    spin_up_year.NH_volcano_size = spin_up_data(row_index,11);
    spin_up_year.SH_volcano_size = spin_up_data(row_index,12);
    if the_year == Param.control_start-1
        write_SO2_flag = 1;
        save_data_flag = 1;
           [Data, emissions] = IAGPperformMPC(the_year, base_directory, project_name, [], [], spin_up_year, save_data_flag, write_SO2_flag);
 %
    else
        write_SO2_flag = 0;
        save_data_flag = 1;
    [Data, emissions] = IAGPperformMPC(the_year, base_directory, project_name, [],[], spin_up_year, save_data_flag, write_SO2_flag);
    end
    row_index = row_index + 1;
end

% do last call to IAGPperformMPC with flag to save data.

