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

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
