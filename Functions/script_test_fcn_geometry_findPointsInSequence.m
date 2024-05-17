% script_test_fcn_geometry_findPointsInSequence
% Exercises the function: fcn_geometry_findPointsInSequence

% Revision history:
% 2023_12_29
% -- wrote the code

%% Set up the workspace

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

assert(length(sequence_indicies) > 1);
assert(length(sequence_indicies(1,:)) == 1);

%% Example - 2 - BASIC call
fig_num = 1;

input_distances = [-1 0 3 4 5 6 7 8.5 9 10 11.5 13 14 15 16 17 18 19 20 21 22]';
base_point_index = 6;
station_tolerance = 2;

sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);

assert(length(sequence_indicies) > 1);
assert(length(sequence_indicies(1,:)) == 1);

%% Example 3
input_distances = [
    0
    0.188896643435862
    0.281186645032550
    0.472304093692704
    0.571973834533390
    0.768818358699212
    0.862778494044884
    0.961593083491137
    1.061653207994030
    1.161395025569664
    10.339478530050306
    10.528375173486166
    10.620665175082857
    10.811782623743010
    10.911452364583695
    11.108296888749514
    11.202257024095186
    11.301071613541442
    11.401131738044334
    11.500873555619966];

fig_num = 2;

base_point_index = 10;
station_tolerance = 0.1;

sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);

assert(length(sequence_indicies) > 1);
assert(length(sequence_indicies(1,:)) == 1);

%% Example - 3 - fast implementation mode
fig_num = 333;

input_distances = [
    0
    0.188896643435862
    0.281186645032550
    0.472304093692704
    0.571973834533390
    0.768818358699212
    0.862778494044884
    0.961593083491137
    1.061653207994030
    1.161395025569664
    10.339478530050306
    10.528375173486166
    10.620665175082857
    10.811782623743010
    10.911452364583695
    11.108296888749514
    11.202257024095186
    11.301071613541442
    11.401131738044334
    11.500873555619966];

base_point_index = 10;
station_tolerance = 0.1;

% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 10; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);
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

%% Example - 2 - fast implementation mode

input_distances = [-1 0 3 6 7 8.5 9 10 11.5 13 14 15 16 19 22]';
base_point_index = 6;
station_tolerance = 2;

% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 1000; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_findPointsInSequence:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Fail cases follow
if 1==0
    %% Fails because the base_point_index is out of range
    input_distances = [
        0
        0.1000
        0.2000
        0.4000
        0.6000
        0.8000
        1.2000
        1.4000
        1.5000
        1.6000
        1.7000
        1.8000
        1.9000
        2.0000
        2.1000
        2.2000
        2.3000
        2.4000
        2.5000
        2.6000
        2.7000
        2.8000
        2.9000
        3.0000
        3.4000
        3.5000
        3.6000
        3.7000
        3.8000
        3.9000
        4.0000
        4.3000
        4.4000
        4.5000
        4.8000
        4.9000
        5.2000
        5.5000
        5.7000
        ];

    fig_num = 4444;

    base_point_index = 91;
    station_tolerance = 1;

    sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance);
end