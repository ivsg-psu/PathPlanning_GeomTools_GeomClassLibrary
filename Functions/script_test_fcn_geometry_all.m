% This is a wrapper script to run all the test scripts in the geometry
% class library for the purpose of evaluating every assertion test in these
% files

% Clear all variables
clear all; close all; clc;

% Find the directory to test
root_directory_name = pwd;

% Does the directory "Functions" exist?
functions_folder_name = fullfile(root_directory_name,'Functions');
if ~exist(functions_folder_name,'dir')
    error('Unable to find the Functions directory')
end

% Make a list of all possible files that start with script_test_fcn
scripts_search_name = fullfile(root_directory_name,'Functions','script_test_fcn_*.m');
all_scripts = dir(scripts_search_name);

% Open a diary to track results
diary 'script_test_fcn_geometry_all_stdout.txt';

% Loop through each file
for i_script = 1:15%length(all_scripts)

    % Get the full name
    file_name_extended = all_scripts(i_script).name;

    % Remove the ".m" extension
    test_script_file_name = erase(file_name_extended,'.m');

    % Check that the file name exists for the script
    if ~strcmp(mfilename,test_script_file_name)        
        suite = testsuite(test_script_file_name);
        results = run(suite);
    end
end
