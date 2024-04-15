%% script_test_all_functions.m
% This is a wrapper script to run all the test scripts in the 
% library for the purpose of evaluating every assertion test in these
% files.
%
% NOTE: to view the output file with formatting, use the "type" command.
% For example:
% type('script_test_fcn_geometry_all_stdout.txt')

clearvars; 
close all; 
clc;
all_scripts = dir(cat(2,'.',filesep,'Functions',filesep,'script_test_fcn_*.m'));
diary 'script_test_fcn_geometry_all_stdout.txt';
for i_script = 1:length(all_scripts)
    file_name_extended = all_scripts(i_script).name;
    file_name = erase(file_name_extended,'.m');
    if ~strcmp(mfilename,file_name)
        file_name_trunc = erase(file_name,'script_');
        fcn_DebugTools_cprintf('*blue',' ');
        fcn_DebugTools_cprintf('*blue','Testing script: %.0d of %.0d, %s\n\n',i_script,length(all_scripts),file_name_trunc);
        suite = testsuite(file_name);
        results = run(suite);
    end
end
diary off