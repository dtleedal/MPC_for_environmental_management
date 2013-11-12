function IAGP_make_SO2_output_file(...
    base_directory,...
    project_name,...
    year,...
    all_SO2_emissions)
% Display and print to file the results
try
    fid = fopen([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'output_file_' num2str(year)],'a+');
    fprintf(fid, '%i\n', year);
    disp('-------------------------------')
    fprintf('%10s%20s\n','Week no.','SO2 emissions (Tg)');
    disp('-------------------------------')
    for week_number = 1:length(all_SO2_emissions)
        fprintf('%7i%15.4f\n', week_number, all_SO2_emissions(week_number));
        fprintf(fid,'%i,%.6f\n', week_number, all_SO2_emissions(week_number));
    end
    disp('-------------------------------')
    fclose(fid);
catch
    error('Could not write the output file')
end


end

