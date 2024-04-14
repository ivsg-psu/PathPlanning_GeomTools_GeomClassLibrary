%% script_test_fcn_geometry_fillArcSequenceTestPoints
% Exercises the function: fcn_geometry_fillArcSequenceTestPoints
% Revision history:
% 2024_03_31 - S. Brennan
% -- wrote the code
% Revision history:
% 2024_04_14 - S. Brennan
% -- added assertions

close all;


%% Test 1: a basic test with 4 points, producing 2 arcs that are similar
fig_num = 1;
figure(fig_num);
clf;

arc_pattern = [1/20, 30; -1/5 9; 0 20; 1/10 50; 0 40; -1/25 39; 0 13; -2/38 38/2*pi];

M = 10;
sigma = 0.02;

[test_points, circleCenters, trueStartPointsOfArcs, arcStartIndicies, namedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, fig_num);

% Check sizes
assert(length(test_points(:,1))>1);
assert(length(test_points(1,:))==2);
assert(length(circleCenters(:,1))==length(arc_pattern(:,1)));
assert(length(circleCenters(1,:))==2);
assert(length(arcStartIndicies(:,1))==length(arc_pattern(:,1)));
assert(length(arcStartIndicies(1,:))==1);
assert(length(namedCurveTypes)==length(arc_pattern(:,1)));
assert(iscell(namedCurveTypes));
assert(length(trueParameters)==length(arc_pattern(:,1)));
assert(iscell(trueParameters));

%% Speed test

% seed_points = [2 3; 4 5; 6 3; 1 1; 2 -1; 0 -3];
% %seed_points = [2 3; 4 5; 6 3; 1 1];
% M = 10;
% sigma = 0.02;
% 
% % Perform the calculation in slow mode
% REPS = 100; minTimeSlow = Inf; 
% tic;
% for i=1:REPS
%     tstart = tic;
%     test_points = fcn_geometry_fillArcSequenceTestPoints(seed_points, M, sigma, []);
%     telapsed = toc(tstart);
%     minTimeSlow = min(telapsed,minTimeSlow);
% end
% averageTimeSlow = toc/REPS;
% 
% % Perform the operation in fast mode
% minTimeFast = Inf;
% tic;
% for i=1:REPS
%     tstart = tic;
%     test_points = fcn_geometry_fillArcSequenceTestPoints(seed_points, M, sigma, -1);
%     telapsed = toc(tstart);
%     minTimeFast = min(telapsed,minTimeFast);
% end
% averageTimeFast = toc/REPS;
% 
% fprintf(1,'\n\nComparison of fcn_geometry_fillArcSequenceTestPoints without speed setting (slow) and with speed setting (fast):\n');
% fprintf(1,'N repetitions: %.0d\n',REPS);
% fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
% fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
% fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
% fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
% fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
% fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end
