function SO2_sum = read_SO2_file( base_directory,...
    project_name,...
    year)

fid = fopen([base_directory filesep project_name filesep 'inputs_and_outputs' filesep 'output_file_' num2str(year)], 'r');
tline = fgetl(fid);
if (str2num(tline) ~= year)
    error('The year in the function call and the year in the file do not coincide')
end
SO2_sum = 0;
while ~feof(fid)
    tline = fgetl(fid);
    
    C = regexp(tline,',','split');
    SO2_sum = SO2_sum + str2num(C{2});
    
end

