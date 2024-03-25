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

%% Test 3: a basic test with 3 points: 180 CCW horizontal
fig_num = 3;


points1 = [0 0];
points2 = [2 -2];
points3 = [4 0];

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,pi));
assert(isequal(circle_centers,[2 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,pi));


%% Test 4: a basic test with 3 points: 180 CW horizontal
fig_num = 4;

points1 = [0 0];
points2 = [2 2];
points3 = [4 0];

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,-pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,-pi));
assert(isequal(circle_centers,[2 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,pi));

%% Test 5: a basic test with 3 points: in same sector
fig_num = 5;

Radius = 2;
points = Radius*[[cos(0)    sin(0)]; [cos(pi/4) sin(pi/4)]; [cos(pi/2) sin(pi/2)]];

points1 = points(1,:);
points2 = points(2,:);
points3 = points(3,:);

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,pi/4));
assert(isequal(arc_angle_in_radians_1_to_3,pi/2));
assert(isequal(circle_centers,[0 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,0));

%% Test 6: a basic test with 3 points: in same sector
fig_num = 6;

Radius = 2;
points = Radius*[[cos(0)    sin(0)]; [cos(pi/4) sin(pi/4)]; [cos(pi/2) sin(pi/2)]];

points1 = points(1,:);
points2 = points(3,:);
points3 = points(2,:);

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,-3*pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,-7*pi/4));
assert(isequal(circle_centers,[0 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,0));

%% Test of fast mode

Radius = 2;
points = Radius*[[cos(0)    sin(0)]; [cos(pi/4) sin(pi/4)]; [cos(pi/2) sin(pi/2)]];

points1 = points(1,:);
points2 = points(3,:);
points3 = points(2,:);

% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
        fcn_geometry_arcAngleFrom3Points(points1, points2, points3, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
        fcn_geometry_arcAngleFrom3Points(points1, points2, points3, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_arcAngleFrom3Points:\n');
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
