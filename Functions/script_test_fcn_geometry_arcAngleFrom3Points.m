% script_test_fcn_geometry_arcAngleFrom3Points
% Exercises the function: fcn_geometry_arcAngleFrom3Points
% Revision history:
% 2023_12_17
% -- wrote the code

close all;
clc;


%% Test 1: a basic test with 3 points: 180 CCW vertical
fig_num = 1;


points1 = [0 0];
points2 = [2 2];
points3 = [0 4];

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,pi));
assert(isequal(circle_centers,[0 2]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,-pi/2));

%% Test 2: a basic test with 3 points: 180 CW vertical
fig_num = 2;

% URHERE

points1 = [0 0];
points2 = [-2 2];
points3 = [0 4];

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,-pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,-pi));
assert(isequal(circle_centers,[0 2]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,-pi/2));

%% Test 3: a basic test with 3 points: 180 CCW vertical
fig_num = 1;


points1 = [0 0];
points2 = [2 2];
points3 = [0 4];

test_points = fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);


%% Test 2: a basic test with 3 points: 180 CW vertical
fig_num = 2;


points1 = [0 0];
points2 = [-2 2];
points3 = [0 4];

test_points = fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);



%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
