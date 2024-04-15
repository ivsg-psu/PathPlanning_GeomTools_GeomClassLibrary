% script_test_fcn_geometry_flagPointsCloserThanLineSegment
% Exercises the function: fcn_geometry_flagPointsCloserToOriginThanLineSegment
% Revision history:
% 2020_06_25 - wrote the code
%
close all;



%% Test 1: a basic test (FAIL)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 -1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 1: a basic test (FAIL - because at origin)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 1: a basic test (PASS)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 1.5];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,1));

%% Test of fast implementation mode 

segment_points = [2 3; 4 5];
test_points = [1 1.5];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Test 1: a basic test (FAIL ABOVE LINE)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [1 4];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 1: a basic test (FAIL ON LINE)
fig_num = 1;
segment_points = [2 3; 4 5];
test_points = [0 1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 2: a basic test vertical line (FAIL)
fig_num = 21;
segment_points = [2 3; 2 5];

test_points = [-1 -1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 2: a basic test vertical line (FAIL - because at origin)
fig_num = 22;
segment_points = [2 3; 2 5];

test_points = [0 0];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 2: a basic test vertical line (PASS)
fig_num = 23;
segment_points = [2 3; 2 5];

test_points = [1 1.5];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,1));

%% Test 2: a basic test vertical line (FAIL ABOVE LINE)
fig_num = 24;
segment_points = [2 3; 2 5];

test_points = [4 4];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 2: a basic test vertical line (FAIL ON LINE)
fig_num = 25;
segment_points = [2 3; 2 5];

test_points = [2 1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 3: a basic test vertical line (FAIL)
fig_num = 31;
segment_points = [-2 3; -2 5];

test_points = [1 -1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 3: a basic test vertical line (FAIL - because at origin)
fig_num = 32;
segment_points = [-2 3; -2 5];

test_points = [0 0];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 3: a basic test vertical line (PASS)
fig_num = 33;
segment_points = [-2 3; -2 5];

test_points = [-1 1.5];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,1));

%% Test 3: a basic test vertical line (FAIL ABOVE LINE)
fig_num = 34;
segment_points = [-2 3; -2 5];

test_points = [-4 4];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 3: a basic test vertical line (FAIL ON LINE)
fig_num = 35;
segment_points = [-2 3; -2 5];

test_points = [-2 1];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

assert(isequal(point_flags,0));

%% Test 4: a basic random test
fig_num = 4;
segment_points = [2 3; 4 5];
test_points = 5*rand(10,2);
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

%% Test 4: a basic random test
fig_num = 3;
segment_points = [1 5; 4 0];
test_points = 5*rand(20,2);
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

%% Test 2: many points randomly determined
fig_num = 21;

segment_points = [0 4; 1 1];
slope = -3;
intercept = 4;

Npoints = 50;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.9*randn(Npoints,1);
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);

%% Test 3: many vertical points
fig_num = 21 + 1;

segment_points = [1 -2; 1  5];
slope = inf;
intercept = inf;

Npoints = 50;
x_data = ones(Npoints,1) + 0.9*randn(Npoints,1);
y_data = linspace(-2,5,Npoints)';
test_points = [x_data,y_data];

[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
% fprintf(1,'Point flags are:\n');
% fprintf(1,'\t%.2f\n',point_flags);


