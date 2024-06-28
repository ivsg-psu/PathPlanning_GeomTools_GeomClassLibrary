%% script_test_fcn_geometry_maxDistanceBetweenPoints
% Exercises the function: fcn_geometry_maxDistanceBetweenPoints

% 2024_04_14 - S. Brennan
% -- wrote the code

close all;

%% Basic example
fig_num = 1;
figure(fig_num);
clf;

test_points_XY = [-1 -1; 2 3; 1 1];
max_distance = fcn_geometry_maxDistanceBetweenPoints(test_points_XY, fig_num);
assert(isequal(round(max_distance,4),round(5,4)));










