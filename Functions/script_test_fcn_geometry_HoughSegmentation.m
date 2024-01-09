% script_test_fcn_geometry_HoughSegmentation
% Exercises the function: fcn_geometry_HoughSegmentation

% Revision history:
% 2023_12_15
% -- wrote the code

close all;
clc;



seed_points = [2 3; 4 5; 7 0; 9 5; 10 20; 13 14];
M = 10; % 10 points per meter
sigma = 0.;
rng(3423)

line_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers to corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

line_segment_test_points = fcn_geometry_corruptPointsWithOutliers(line_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption));

%% Basic example 1: find 5 lines within noisy data
fig_num = 1;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = inf; % Units are meters
threshold_max_points = 10;
input_points = line_segment_test_points;
domains = fcn_geometry_HoughSegmentation(line_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example 2: find 5 line segments within same data
% segments are created by imposing constraints on separation

% Call the segmentation function
fig_num = 2;
transverse_tolerance = 0.02; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = line_segment_test_points;
domains = fcn_geometry_HoughSegmentation(line_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example 3: find segments within a chevron
M = 10; % 40 points per meter

rng(234)
sigma = 0.02;

line_segment_test_points = [];

seed_points = [0 0; 10 0];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points]; %#ok<*NASGU>

seed_points = [0 0; 10 5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

seed_points = [2 0; 3 1.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

seed_points = [4 0; 5 2.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

seed_points = [6 0; 7 3.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

seed_points = [8 0; 9 4.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

seed_points = [10 0; 10 5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
line_segment_test_points = [line_segment_test_points; subtest_points];

% Add outliers to corrupt the results
outliers = [10*rand(100,1) 5*rand(100,1)];
line_segment_test_points = [line_segment_test_points; outliers];


% Call the segmentation function
fig_num = 3;
transverse_tolerance = 0.02; % Units are meters
station_tolerance = 0.6; % Units are meters
threshold_max_points = 10;
input_points = line_segment_test_points;
domains = fcn_geometry_HoughSegmentation(line_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


%% Test of fast mode
% Perform the calculation in slow mode
REPS = 1; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    domains = fcn_geometry_HoughSegmentation(line_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    domains = fcn_geometry_HoughSegmentation(line_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_HoughSegmentation:\n');
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
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

