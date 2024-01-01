% script_test_fcn_geometry_findPointsInSequence
% Exercises the function: fcn_geometry_findPointsInSequence

% Revision history:
% 2023_12_29
% -- wrote the code

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

%% Example - 1 - BASIC call
fig_num = 1;

input_distances = [-1 0 3 6 7 8.5 9 10 11.5 13 14 15 16 19 22]';
base_point_index = 6;
station_tolerance = 2;

sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);

%% Example - 2 - fast implementation mode
fig_num = 1;

inputPoints = test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;

% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [fittedParameters, agreementIndices] = fcn_geometry_findPointsInSequence(inputPoints, transverse_tolerance, station_tolerance, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 10; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [fittedParameters, agreementIndices] = fcn_geometry_findPointsInSequence(inputPoints, transverse_tolerance, station_tolerance, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_findPointsInSequence:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);


%% Example 3
input_distances = [
    0
    1.1069
    0.3847
    1.5977
    4.6230
    0.6998
    2.3057
    0.1035
    5.1142
    2.1066
    0.9996
    2.7229
    0.8983
    4.8260
    1.8245
    0.8054
    2.2122
    1.3189
    0.2847
    0.4960
    -0.0956
    0.2146
    -0.3027
    1.4045
    1.6931
    3.9222
    -0.1815
    -0.5896];

fig_num = 3;

base_point_index = 1;
station_tolerance = 0.4;

sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);
