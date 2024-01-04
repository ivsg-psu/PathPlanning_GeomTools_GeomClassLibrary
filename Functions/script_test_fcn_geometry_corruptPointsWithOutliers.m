% script_test_fcn_geometry_corruptPointsWithOutliers
% Exercises the function: fcn_geometry_corruptPointsWithOutliers
% Revision history:
% 2023_12_05
% -- wrote the code

close all;
clc;


%% Test 1: a basic test with just 2 points, few outliers
fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the results
probability_of_corruption = [];
magnitude_of_corruption = [];

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));



%% Test 2: a basic test with just 3 points
fig_num = 2;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 7 0];
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the results
probability_of_corruption = [];
magnitude_of_corruption = [];

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


%% Test 3: a basic test with large corruption
fig_num = 3;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 7 0];
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

%% Test of fast mode

% Perform the calculation in slow mode
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_corruptPointsWithOutliers:\n');
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
