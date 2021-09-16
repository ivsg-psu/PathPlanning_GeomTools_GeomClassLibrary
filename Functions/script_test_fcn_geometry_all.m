% This is a wrapper script to run all the test scripts in the geometry
% class library for the purpose of evaluating every assertion test in these
% files
clear all; close all; clc;
all_scripts = dir('script_test_fcn_*');
suites = [];

for i_script = 1:2%length(all_scripts)
    file_name_extended = all_scripts(i_script).name;
    file_name = erase(file_name_extended,'.m');
    if ~strcmp(mfilename,file_name) && ~strcmp(file_name(end-3:end),'.asv')
        file_name_trunc = erase(file_name,'script_');
        suite = testsuite(file_name);
        suites(end+1) = suite;
    end
end
results = run(suites)
