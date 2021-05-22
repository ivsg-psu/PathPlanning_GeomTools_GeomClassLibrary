% script_test_fcn_geometry_flagPointsCloserThanLineSegment
% Exercises the function: fcn_geometry_flagPointsCloserToOriginThanLineSegment
% Revision history:
% 2020_06_25 - wrote the code
%
close all;
clc;



%% Test 1: a basic test
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: many points randomly determined
fig_num = fig_num + 1;

segment_points = [0 4; 1 1];
slope = -3;
intercept = 4;

Npoints = 50;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: many vertical points
fig_num = fig_num + 1;

segment_points = [1 -2; 1  5];
slope = inf;
intercept = inf;

Npoints = 50;
x_data = ones(Npoints,1) + 0.2*randn(Npoints,1);
y_data = linspace(-2,5,Npoints)';
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);


