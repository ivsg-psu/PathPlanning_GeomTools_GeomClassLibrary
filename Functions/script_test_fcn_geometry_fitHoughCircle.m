% script_test_fcn_geometry_circleHoughFit
% Exercises the function: fcn_geometry_circleHoughFit

% Revision history:
% 2023_12_15
% -- wrote the code
% 2024_01_04
% -- fixed the argument inputs

%% Set up the workspace

clc
close all

%% Examples for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ______                           _
% |  ____|                         | |
% | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                            | |
%                            |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Fill test data - 1 segment
fig_num = 23;
% figure(fig_num);
% clf;
% hold on;
% axis equal
% grid on;

seed_points = [2 3; 4 5; 6 3]; %; 1 1];
M = 5;
sigma = 0.02;

test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

%% Example - 1 - BASIC call
fig_num = 1;

inputPoints = test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, fig_num);

%% Example - 2 - fast implementation mode
fig_num = 2222;

inputPoints = test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;

% Perform the calculation in slow mode
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 10; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitHoughCircle:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

