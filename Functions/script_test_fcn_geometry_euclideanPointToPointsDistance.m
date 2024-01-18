% script_test_fcn_geometry_euclideanPointToPointsDistance
% Tests fcn_geometry_euclideanPointToPointsDistance

% Revision History:
% 2024_01_18 - Aneesh Batchu

%% BASIC example - two points in 2D
fig_num = 2;

pt1 = [1 1; 0 0];
pt2 = [2 3; 4 0] ;
dist=fcn_geometry_euclideanPointToPointsDistance(pt1,pt2,fig_num);
assert(isequal(round(dist,4), [2.2361; 4]));

%% BASIC example - many points in 2D
fig_num = 3;

pt1 = rand(5,2);
pt2 = rand(5,2);
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

% BASIC example - single points in 3D
fig_num = 31;

pt1 = [1 1 0];
pt2 = [2 3 2];
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);
assert(isequal(round(dist,4), [3]));

%% BASIC example - two points in 3D
fig_num = 32;

pt1 = [1 1 0; 0 0 1];
pt2 = [2 3 4; 4 0 2] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);
assert(isequal(round(dist,4), [4.5826; 4.1231]));

%% BASIC example - multiple points in 3D
fig_num = 33;

pt1 = [-1 1 0; 0 0 1; -3 -2 -4];
pt2 = [2 3 4; 4 0 2; -5 3 -2] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);
assert(isequal(round(dist,4), [5.3852; 4.1231; 5.7446]));

%% BASIC example - multiple points in 3D
fig_num = 34;

pt1 = rand(5,3);
pt2 = rand(5,3);
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% Test of fast implementation mode 

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    dist = fcn_geometry_euclideanPointsToPointsDistance(pt1, pt2, (fig_num));
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
    dist = fcn_geometry_euclideanPointsToPointsDistance(pt1, pt2, (fig_num));
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

