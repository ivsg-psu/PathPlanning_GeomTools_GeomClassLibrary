%% script_test_fcn_geometry_pointsNearPoint
% This is a script to exercise the function:
% fcn_geometry_pointsNearPoint.m
% This function was written on 2024_08_26 by S. Brennan, sbrennan@psu.edu


%% test 1 - Simple example
fig_num = 1;
figure(fig_num);
clf;

% Fill some data
anchorPoint = [0 0];
pointsToSearch = [linspace(0,5,6)' zeros(6,1)];
searchRadius = 2;

% Test the function
nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_pointsNearPoint',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing simple example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 3)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 3]'));



%% test 2 - Simple example - no plotting
fig_num = 2;
figure(fig_num);
close(fig_num);

% Fill some data
anchorPoint = [0 0];
pointsToSearch = [linspace(0,5,6)' zeros(6,1)];
searchRadius = 2;

% Test the function
nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, ([]));


% Was a figure created?
assert(all(~ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 3)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 3]'));


%% test 3 - empty result example
fig_num = 3;
figure(fig_num);
clf;

% Fill some data
anchorPoint = [0 0];
pointsToSearch = [linspace(4,5,3)' zeros(3,1)];
searchRadius = 2;

% Test the function
nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_pointsNearPoint',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing empty result example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(isempty(nearbyIndicies));


%% test 4 - 3D points example
fig_num = 4;
figure(fig_num);
clf;

% Fill some data
anchorPoint = [0 0 0];
pointsToSearch = [linspace(0,5,6)' zeros(6,2)];
searchRadius = 2;

% Test the function
nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_pointsNearPoint',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing 3D points example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 3)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 3]'));

%% test 5 - Lots of points
fig_num = 5;
figure(fig_num);
clf;

% Fill some data
anchorPoint = [2 3];
pointsToSearch = 10*rand(50,2) - ones(50,2)*5;
searchRadius = 2;

% Test the function
nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_pointsNearPoint',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing many points example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
%assert(length(nearbyIndicies(:,1))== 3)

% % Does the data contain the correct values?
% assert(isequal(nearbyIndicies,[1 2 3]'));

%% Speed test

% Fill some data
anchorPoint = [0 0 0];
pointsToSearch = [linspace(0,5,6)' zeros(6,2)];
searchRadius = 2;


% Test the function
fig_num=[];
REPS=5; 
minTimeSlow=Inf;
maxTimeSlow=-Inf;
tic;

% Slow mode calculation - code copied from plotVehicleXYZ
for i=1:REPS
    tstart=tic;
    nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
    telapsed=toc(tstart);
    minTimeSlow=min(telapsed,minTimeSlow);
    maxTimeSlow=max(telapsed,maxTimeSlow);
end
averageTimeSlow=toc/REPS;
% Slow mode END

% Fast Mode Calculation
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;
% Fast mode END

% Display Console Comparison
if 1==1
    fprintf(1,'\n\nComparison of fcn_geometry_pointsNearPoint without speed setting (slow) and with speed setting (fast):\n');
    fprintf(1,'N repetitions: %.0d\n',REPS);
    fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
    fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
    fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
    fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
    fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
    fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',maxTimeSlow/minTimeFast);
end
%Assertion on averageTime NOTE: Due to the variance, there is a chance that
%the assertion will fail.
assert(averageTimeFast<4*averageTimeSlow);

