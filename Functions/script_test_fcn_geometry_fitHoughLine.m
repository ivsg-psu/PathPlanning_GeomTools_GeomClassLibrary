% script_test_fcn_geometry_fitHoughLine
% Exercises the function: fcn_geometry_fitHoughLine
% Revision history:
% 2023_12_14
% -- wrote the code

close all;
clc;


%% Test 1: a basic test
fig_num = 1;
transverse_tolerance = 0.2;

[fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine([1 0; 1 1], transverse_tolerance, station_tolerance, fig_num);


%% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

%% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));

%% Test 2: a basic test
fig_num = 2;
figure(fig_num);
clf;

transverse_tolerance = 0.2;

[fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance,  fig_num);

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitHoughLine(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end