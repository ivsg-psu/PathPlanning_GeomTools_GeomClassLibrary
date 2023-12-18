% script_test_fcn_geometry_shufflePointOrdering
% Exercises the function: fcn_geometry_shufflePointOrdering
% Revision history:
% 2023_12_17
% -- wrote the code

close all;
clc;


%% Test 1: a basic test 
fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 3 2];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

shuffled_points = fcn_geometry_shufflePointOrdering(test_points, (fig_num));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
