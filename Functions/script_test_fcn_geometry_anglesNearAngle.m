%% script_test_fcn_geometry_anglesNearAngle
% This is a script to exercise the function:
% fcn_geometry_anglesNearAngle.m
% This function was written on 2024_08_26 by S. Brennan, sbrennan@psu.edu


%% test 1 - Simple example
fig_num = 1;
figure(fig_num);
clf;

% Fill some data
anchorAngle = 0;
anglesToSearch = linspace(0,360,37)'*pi/180;
angleRange = 15*pi/180;

% Test the function
nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_anglesNearAngle',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing simple example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 4)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 36 37]'));



%% test 2 - Simple example - no plotting
fig_num = 2;
figure(fig_num);
close(fig_num);

% Fill some data
anchorAngle = 0;
anglesToSearch = linspace(0,360,37)'*pi/180;
angleRange = 15*pi/180;

% Test the function
nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, ([]));


% Was a figure created?
assert(all(~ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 4)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 36 37]'));


%% test 3 - empty result example
fig_num = 3;
figure(fig_num);
clf;

% Fill some data
anchorAngle = 90;
anglesToSearch = linspace(0,40,6)'*pi/180;
angleRange = 15*pi/180;

% Test the function
nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_anglesNearAngle',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing empty result example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(isempty(nearbyIndicies));


%% test 4 - rollover angles example
fig_num = 4;
figure(fig_num);
clf;

% Fill some data
anchorAngle = 180*pi/180;
anglesToSearch = linspace(-180,180,37)'*pi/180;
angleRange = 15*pi/180;

% Test the function
nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num));
title(sprintf('Example %.0d: showing fcn_geometry_anglesNearAngle',fig_num), 'Interpreter','none','FontSize',12);
subtitle('Showing rollover angles example');

% Was a figure created?
assert(all(ishandle(fig_num)));

% Does the data have right size?
assert(~isempty(nearbyIndicies));
assert(1==length(nearbyIndicies(1,:)));
assert(length(nearbyIndicies(:,1))== 4)

% Does the data contain the correct values?
assert(isequal(nearbyIndicies,[1 2 36 37]'));

%% Speed test

% Load the data
anchorAngle = 0;
anglesToSearch = linspace(0,360,37)'*pi/180;
angleRange = 15*pi/180;


% Test the function
fig_num=[];
REPS=5; 
minTimeSlow=Inf;
maxTimeSlow=-Inf;
tic;

% Slow mode calculation - code copied from plotVehicleXYZ
for i=1:REPS
    tstart=tic;
    nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num));
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
    nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;
% Fast mode END

% Display Console Comparison
if 1==1
    fprintf(1,'\n\nComparison of fcn_geometry_anglesNearAngle without speed setting (slow) and with speed setting (fast):\n');
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

