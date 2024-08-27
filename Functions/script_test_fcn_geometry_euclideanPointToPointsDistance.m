% script_test_fcn_geometry_euclideanPointToPointsDistance
% Tests fcn_geometry_euclideanPointToPointsDistance

% Revision History:
% 2024_01_18 - Aneesh Batchu
% -- added basic assertions

%% BASIC example - two points in 2D
fig_num = 201;
figure(fig_num);
clf;

pt1 = [1 1];
pt2 = [1 3; 4 1] ;
dist=fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, fig_num);
title(sprintf('Example %.0d: showing fcn_geometry_euclideanPointToPointsDistance',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing 2D points');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(1==length(dist(1,:)));
assert(length(dist(:,1))== length(pt2(:,1)))

% Does the data contain the correct values?
assert(isequal(round(dist,4), [2 3]'));

%% BASIC example - two points in 2D, no plotting
fig_num = 202;
figure(fig_num);
close(fig_num)

pt1 = [1 1];
pt2 = [1 3; 4 1] ;
dist=fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, []);

% Was a figure created?
assert(all(~ishandle(fig_num)));

% Does the data have right size?
assert(1==length(dist(1,:)));
assert(length(dist(:,1))== length(pt2(:,1)))

% Does the data contain the correct values?
assert(isequal(round(dist,4), [2 3]'));

%% BASIC example - many points in 2D
fig_num = 203;
figure(fig_num);
clf;

pt1 = [0 0];
angles = 2*pi*rand(6,1);
pt2 = [cos(angles) sin(angles)];
dist=fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, fig_num);
title(sprintf('Example %.0d: showing fcn_geometry_euclideanPointToPointsDistance',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing many random 2D points');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(1==length(dist(1,:)));
assert(length(dist(:,1))== length(pt2(:,1)))

% Does the data contain the correct values?
assert(isequal(round(dist,4), ones(length(angles),1)));


%% BASIC example - two points in 3D
fig_num = 301;
figure(fig_num);
clf;

pt1 = [1 1 0];
pt2 = [1 1 4; 4 1 0] ;
dist=fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, fig_num);
title(sprintf('Example %.0d: showing fcn_geometry_euclideanPointToPointsDistance',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing 3D points');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(1==length(dist(1,:)));
assert(length(dist(:,1))== length(pt2(:,1)))

% Does the data contain the correct values?
assert(isequal(round(dist,4), [4 3]'));


%% Test of fast implementation mode 

pt1 = rand(1,3);
pt2 = rand(5,3);

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    dist = fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, (fig_num));
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
    dist = fcn_geometry_euclideanPointToPointsDistance(pt1, pt2, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_euclideanPointToPointsDistance:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

