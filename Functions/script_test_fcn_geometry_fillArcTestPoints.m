% script_test_fcn_geometry_fillArcTestPoints
% Exercises the function: fcn_geometry_fillArcTestPoints
% Revision history:
% 2023_12_17
% -- wrote the code

close all;
clc;


%% Test 1: a basic test with 4 points, producing 2 arcs that are similar
fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 6 3; 1 1];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
