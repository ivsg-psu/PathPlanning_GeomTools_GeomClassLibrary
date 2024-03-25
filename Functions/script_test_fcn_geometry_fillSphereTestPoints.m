% script_test_fcn_geometry_fillSphereTestPoints
% Exercises the function: fcn_geometry_fillSphereTestPoints
% Revision history:
% 2024_01_23 - S. Brennan
% -- wrote the code

close all;
clc;


%% Test 1: a basic test with 200 points
fig_num = 1;
figure(fig_num);
clf;

N_points = 200;
sphere_center = [3 5 0];
sphere_radius = 2;
sigma = 0.02;

test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num));
assert(length(test_points(:,1))==N_points);


%% Test 1: a basic test with 100 points
fig_num = 1;

N_points = 100;
sphere_center = [6 2 0];
sphere_radius = 1;
sigma = 0.02;

test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num));
assert(length(test_points(:,1))==N_points);


%% Speed test

N_points = 100;
sphere_center = [6 2 0];
sphere_radius = 1;
sigma = 0.02;

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fillSphereTestPoints without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
