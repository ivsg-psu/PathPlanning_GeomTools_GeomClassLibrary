% script_test_fcn_geometry_arcDirectionFrom3Points - this is is a script written
% to test the function: fcn_geometry_arcDirectionFrom3Points.m.
%

% Revision history:
% 2023_12_19 - sbrennan@psu.edu
% -- original write of the code
close all;

%% BASIC example 1 - positive
fig_num = 1;
figure(fig_num);clf;

points1 = [0 0];
points2 = [1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,1))

%% BASIC example 2 - negative
fig_num = 2;
figure(fig_num);clf;

points1 = [0 0];
points2 = [-1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,-1))

%% BASIC example 3 - stacked
fig_num = 3;
figure(fig_num);clf;

points1 = [0 0; 3 3;];
points2 = [-1 4; 4 4];
points3 = [0 5; 3 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,[-1; 1]))

%% Test of fast mode

% Perform the calculation in slow mode
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_arcDirectionFrom3Points:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Debugging case
fig_num = 1232;
figure(fig_num);clf;

seed_points = [2 3; 4 5; 6 3; 1 1; -2 -1];

% Find if the arcs are counterclockwise or clockwise
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(seed_points(1:end-2,:), seed_points(2:end-1,:), seed_points(3:end,:), fig_num);
assert(isequal(is_counterClockwise,[-1; -1; 1]))


%% Fail cases follow
if 1==0
    %% FAIL 1: points2 not same length
    points1 = [2 3; 3 4];
    points2 = [4 5];
    points3 = [5 6; 7 8];
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3,fig_num);

end



