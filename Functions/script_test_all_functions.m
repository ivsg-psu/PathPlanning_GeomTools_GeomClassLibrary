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

N_files = length(all_scripts);
testing_times = nan(N_files,1);

diary 'script_test_fcn_geometry_all_stdout.txt';
dbclear if warning
% dbstop if warning

%  badScripts = [3 5 6 9 10 11 17 23 26 27 55 81];
badScriptNames = {
    'script_test_fcn_geometry_HoughRegression';          % 3   <--Aneesh
    'script_test_fcn_geometry_alignArcArc';              % 5   <--Sean
    'script_test_fcn_geometry_alignArcArcC2Optimized';   % 6   <--Sean
    'script_test_fcn_geometry_alignArcsInSequence';      % 9   <--Sean
    'script_test_fcn_geometry_alignGeometriesInSequence';% 10  <--Sean
    'script_test_fcn_geometry_alignLineArc';             % 11  <--Sean
    'script_test_fcn_geometry_boundaryAnalysis';         % 17  <--Aneesh
    'script_test_fcn_geometry_compareCurves';            % 23  <--Sean
    'script_test_fcn_geometry_concatenatePoints';        % 26  <--Aneesh
    'script_test_fcn_geometry_concentricCubesPointDensity'; % 27  <--Aneesh
    'script_test_fcn_geometry_findDrivenPath';           % 55  <--Aneesh
    'script_test_fcn_geometry_fitSequentialArcs';        % 81  <--Sean
    'script_test_fcn_geometry_intersectGeom';            % 91  <--Sean
    'script_test_fcn_geometry_isC1FeasibleArcToArc';     % 93  <--Sean
    'script_test_fcn_geometry_isC2FeasibleArcToArc';     % 94  <--Sean
    'script_test_fcn_geometry_isFeasibleAlignGeomPair';  % 97  <--Sean
    'script_test_fcn_geometry_isFeasibleAlignGeomSeries';% 98  <--Sean
    'script_test_fcn_geometry_isFeasibleGeomSequence';   % 99  <--Sean
    'script_test_fcn_geometry_plotGeometry';             % 108 <--Sean
    'script_test_fcn_geometry_stdInZ';                   % 119 <--Aneesh
    'script_test_fcn_geometry_surfaceAnalysis';          % 120 <--Aneesh
    'script_test_fcn_geometry_vehiclePosition';          % 120 <--Aneesh
    };
badScripts = (1:N_files);

for scriptIndex = 1:length(badScripts)
    i_script = badScripts(scriptIndex);
    file_name_extended = all_scripts(i_script).name;
    file_name = erase(file_name_extended,'.m');
    if ~strcmp(mfilename,file_name)
        %file_name_trunc = erase(file_name,'script_');
        fcn_DebugTools_cprintf('*blue',' ');
        fcn_DebugTools_cprintf('*blue','Testing script: %.0d of %.0d, %s\n\n',i_script,length(all_scripts),file_name);
        % disp('Press any key to continue');
        tstart = tic;
        suite = testsuite(file_name);
        results = run(suite);
        telapsed = toc(tstart);
        testing_times(i_script) = telapsed;
    end
end
diary off

close all;
figure(458908);
plot(testing_times);
grid on;
xlabel('Script test number');
ylabel('Elapsed time to test (sec)');

fprintf(1,'The testing times for each script:\n');
for i_script = 1:N_files
    if ~isnan(testing_times(i_script))
        fprintf(1,'%.0d: \t %.2f seconds for script: \t %s \n',i_script, testing_times(i_script), all_scripts(i_script).name);
    end
end