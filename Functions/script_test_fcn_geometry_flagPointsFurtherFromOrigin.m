% script_test_fcn_geometry_flagPointsFurtherFromOrigin.m
% Exercises the function: fcn_geometry_flagPointsFurtherFromOriginThanLineSegment

% Revision history:
% 2021_05_31
% -- Wrote the code via mods to fcn_geometry_flagPointsCloserToOriginThanLineSegment

close all;
clc;



%% Test 1: a basic test (FAIL)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 -1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 1: a basic test (FAIL - because at origin)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 1: a basic test (FAIL)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 1.5];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 1: a basic test (PASS ABOVE LINE)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 4];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 1: a basic test (FAIL ON LINE)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [0 1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);




%% Test 2: a basic test vertical line (FAIL)
fig_num = 21;
segment_points = [2 3; 2 5];

test_points = [-1 -1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: a basic test vertical line (FAIL - because at origin)
fig_num = 22;
segment_points = [2 3; 2 5];

test_points = [0 0];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: a basic test vertical line (FAIL)
fig_num = 23;
segment_points = [2 3; 2 5];

test_points = [1 1.5];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: a basic test vertical line (PASS ABOVE LINE)
fig_num = 24;
segment_points = [2 3; 2 5];

test_points = [4 4];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: a basic test vertical line (FAIL ON LINE)
fig_num = 25;
segment_points = [2 3; 2 5];

test_points = [2 1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);







%% Test 3: a basic test vertical line (FAIL)
fig_num = 31;
segment_points = [-2 3; -2 5];

test_points = [1 -1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: a basic test vertical line (FAIL - because at origin)
fig_num = 32;
segment_points = [-2 3; -2 5];

test_points = [0 0];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: a basic test vertical line (FAIL)
fig_num = 33;
segment_points = [-2 3; -2 5];

test_points = [-1 1.5];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: a basic test vertical line (PASS ABOVE LINE)
fig_num = 34;
segment_points = [-2 3; -2 5];

test_points = [-4 4];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: a basic test vertical line (FAIL ON LINE)
fig_num = 35;
segment_points = [-2 3; -2 5];

test_points = [-2 1];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);





%% Test 4: a basic random test
fig_num = 41;
segment_points = [2 3; 4 5];
test_points = 5*rand(10,2);
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 4: a basic random test
fig_num = 42;
segment_points = [1 5; 4 0];
test_points = 5*rand(20,2);
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 4: many points randomly determined around line
fig_num = 43;

segment_points = [0 4; 1 1];
slope = -3;
intercept = 4;

Npoints = 50;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.9*randn(Npoints,1);
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: many vertical points around random vertical line
fig_num = 44

segment_points = [1 -2; 1  5];
slope = inf;
intercept = inf;

Npoints = 50;
x_data = ones(Npoints,1) + 0.9*randn(Npoints,1);
y_data = linspace(-2,5,Npoints)';
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);


