% script_test_fcn_geometry_calcUnitVector
% Exercises the function: fcn_geometry_calcUnitVector
% Revision history:
% 2023_12_14
% -- wrote the code

close all;
clc;


%% Test 1: a basic test
fig_num = 1;
input_vectors = [3 3]; 

unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num); %#ok<NASGU>


%% Test 2: many vectors
fig_num = 2;
input_vectors = randn(10,2); 
unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

%% Test 3: fast mode
% Perform the calculation in slow mode
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    unit_vectors = fcn_geometry_calcUnitVector(input_vectors, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 1000; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    unit_vectors = fcn_geometry_calcUnitVector(input_vectors, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_calcUnitVector:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_calcUnitVector(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end