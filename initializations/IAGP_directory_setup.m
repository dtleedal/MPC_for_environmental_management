function IAGP_directory_setup(base_directory,project_name)
% script to generate model parameter sets
% D.T. Leedal
% 27/2/2013
%
% define the base working directory
% NOTE TO LAWRENCE:
% you need to change this to the path where *you* placed MPC_for_Lawrence directory.
% This path will always end with \MPC_for_Lawrence
% If you are already in the MPC_for_Lawrence directory you can type pwd in
% the command window and copy and paste the result (but don't forget to
% surround with apostrophes.
%base_directory = ['C:' filesep 'Documents and Settings' filesep 'leedald' filesep 'My Documents' filesep 'Dropbox' filesep 'IAGP_HadGEM_control' filesep 'AJs_control_model' filesep 'MPC_for_Lawrence'];
% for example on my Mac it would be:
%base_directory = [filesep 'Users' filesep 'davidleedal' filesep 'Dropbox' filesep 'IAGP_HadGEM_control' filesep 'AJs_control_model' filesep 'MPC_for_Lawrence'];

% give the project a name to use as the directory for the project
% best to avoid spaces in the name
%project_name = 'dave_test_project';
disp(['project name: '  project_name])

%% make the project directory structure
disp('Making directories...')
try
    mkdir([base_directory filesep project_name])
    addpath([base_directory filesep project_name])
    % make the experiment state data directory
    mkdir([base_directory filesep project_name filesep 'state_data']);
    addpath([base_directory filesep project_name filesep 'state_data']);
    disp(['Made: ' filesep project_name filesep 'state_data'])
    mkdir([base_directory filesep project_name filesep 'backup']);
    addpath([base_directory filesep project_name filesep 'backup']);
    disp(['Made: ' filesep project_name filesep 'backup'])
    mkdir([base_directory filesep project_name filesep 'inputs_and_outputs']);
    addpath([base_directory filesep project_name filesep 'inputs_and_outputs']);
    disp(['Made: ' filesep project_name filesep 'inputs_and_outputs'])
    mkdir([base_directory filesep project_name filesep 'figures']);
    addpath([base_directory filesep project_name filesep 'figures']);
    disp(['Made: ' filesep project_name filesep 'figures'])
    disp('...done')
catch
    error('Could not write 1 or more of the necessary directories...')
end
