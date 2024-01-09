% script_test_fcn_geometry_fitSlopeInterceptNPoints
% Exercises the function: fcn_geometry_fitSlopeInterceptNPoints
% Revision history:
% 2020_06_25
% -- wrote the code
% 2021_05_24
% -- revised name to fcn_geometry_fitSlopeInterceptNPoints

close all;
clc;
fig_num = 1;

%% Test 1: a basic test
points = [2 3; 4 5];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
assert(isequal(slope,1));
assert(isequal(intercept,1));

%% Test 2: a vertical line
fig_num = fig_num + 1;
points = [2 0; 2 2];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 3: many points randomly determined
fig_num = fig_num + 1;
slope = -3;
intercept = 4;
Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
points = [x_data,y_data];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 4: many vertical points
fig_num = fig_num + 1;
Npoints = 1000;
x_data = 2*ones(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 5: a singular situation (determinant is zero - which gives b = 0)
fig_num = fig_num + 1;
points = [6 4; 3 2];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test of fast mode
% Fill in data to test
slope = -3;
intercept = 4;
Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
points = [x_data,y_data];


% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 1000; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitSlopeInterceptNPoints:\n');
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
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end