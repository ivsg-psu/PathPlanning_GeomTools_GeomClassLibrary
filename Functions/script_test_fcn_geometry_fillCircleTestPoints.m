% script_test_fcn_geometry_fillCircleTestPoints
% Exercises the function: fcn_geometry_fillCircleTestPoints
% Revision history:
% 2023_12_17
% -- wrote the code

close all;


%% Test 1: a basic test with 4 points, producing 2 arcs that are similar
fig_num = 1;
figure(fig_num);
clf;

circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;

test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma, (fig_num));


%% Speed test

circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma, (fig_num));
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
    test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma, (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fillCircleTestPoints without speed setting (slow) and with speed setting (fast):\n');
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
