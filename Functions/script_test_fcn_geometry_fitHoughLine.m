% script_test_fcn_geometry_fitHoughLine
% Exercises the function: fcn_geometry_fitHoughLine
% Revision history:
% 2023_12_14
% -- wrote the code

close all;

%% Fill in some test data
rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);


%% Test 1: a basic test of line segment fitting
fig_num = 1;
transverse_tolerance = 0.2;
station_tolerance = 2;
points_required_for_agreement = [];

domains= fcn_geometry_fitHoughLine([1 0; 1 1], transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);


%% Test 2: a basic test of line segment fitting, noisy points

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

fig_num = 2;
figure(fig_num);
clf;

transverse_tolerance = 0.2;
station_tolerance = 0.4;
points_required_for_agreement = [];

domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);

%% Test 3: a basic test of LINE fitting, noisy points, showing effect of station_tolerance setting

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

fig_num = 3;
figure(fig_num);
clf;

transverse_tolerance = 0.2;
station_tolerance = [];
points_required_for_agreement = [];

domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement,  fig_num);

%% Test 4: a basic test of line segment fitting, noisy points

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

fig_num = 4;
figure(fig_num);
clf;

transverse_tolerance = 0.2;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);

%% Test 5: confirming that it does not plot if figure is off

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

fig_num = 5;

transverse_tolerance = 0.2;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains = fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement);

fcn_geometry_plotFitDomains(domains, fig_num);

%% Test 6: show get empty domain if points for agreement is too high

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

fig_num = 6;
figure(fig_num);
clf;

transverse_tolerance = 0.2;
station_tolerance = 0.4;
points_required_for_agreement = length(test_points(:,1))+1;

domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);


%% Test of fast mode

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

% Set values for testing
transverse_tolerance = 0.2;
station_tolerance = 1;
points_required_for_agreement = [];


% Perform the calculation in slow mode
REPS = 10; 
minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    domains= ...
        fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement,  []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitHoughLine:\n');
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
    [slope,intercept] = fcn_geometry_fitHoughLine(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end