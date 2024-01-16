% script_test_fcn_geometry_findAgreementsOfPointsToCircle
% Exercises the function: fcn_geometry_findAgreementsOfPointsToCircle

% Revision history:
% 2024_01_15 - S. Brennan
% -- wrote the code

close all;
clc;

%% Fill in test data
rng(383);


fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Fill test data for 2 arcs
first_fraction = [0 0.5]; % data from 0 to 50 percent
second_fraction = [0.80 1]; % data from 80 percent to end
N_points = length(onearc_test_points(:,1));

first_fraction_indicies = round(first_fraction*N_points); % find closest indicies
first_fraction_indicies = max([first_fraction_indicies; 1 1],[],1); % Make sure none are below 1
first_fraction_indicies = min([first_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

second_fraction_indicies = round(second_fraction*N_points); % find closest indicies
second_fraction_indicies = max([second_fraction_indicies; 1 1],[],1); % Make sure none are below 1
second_fraction_indicies = min([second_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

twoarc_test_points = ...
    [onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];

corrupted_twoarc_test_points = ...
    [corrupted_onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    corrupted_onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];


% For debugging
figure(33838);
clf;
hold on;
grid on;
grid minor;
axis equal;
plot(corrupted_twoarc_test_points(:,1),corrupted_twoarc_test_points(:,2),'k.');

% 1 outlier arc
seed_points = [6 6; 9 3; 6 0];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 8; % Number of points per meter
sigma = 0.02;

outlieronearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

figure(234);
clf;
hold on;
grid on;
grid minor;
axis equal;

corrupted_outlieronearc_test_points= fcn_geometry_corruptPointsWithOutliers(outlieronearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (234));




%% Test 1: a basic test of circle fitting
fig_num = 1;
figure(fig_num); clf;

points = corrupted_outlieronearc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
transverse_tolerance = 0.2;


agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num));


%% Test of fast mode
% Set inputs to a "normal" mode
base_point_index = 1;
end_point = [9 4];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 1;
station_tolerance = 2;

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num));
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
    agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findAgreementsOfPointsToCircle:\n');
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
    [slope,intercept] = fcn_geometry_findAgreementsOfPointsToCircle(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end