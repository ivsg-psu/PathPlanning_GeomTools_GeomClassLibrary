% script_test_fcn_geometry_fitArcRegressionFromHoughFit
% Exercises the function: fcn_geometry_fitArcRegressionFromHoughFit

% Revision history:
% 2024_01_09 - S. Brennan
% -- wrote the code

close all;
clc;


%% Filling test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];
corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];

start_vector = arc_seed_points(1,:)-arc_true_circleCenter;
arc_true_start_angle_in_radians = atan2(start_vector(2),start_vector(1));
end_vector = arc_seed_points(end,:)-arc_true_circleCenter;
arc_true_end_angle_in_radians = atan2(end_vector(2),end_vector(1));

% Fill circle data
% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));
circle_true_parameters = [circle_center, circle_radius, 0, 2*pi, 1];

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));



%% Fit the onarc_test_points
fig_num = 234;
figure(fig_num); clf;

transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = 20;
flag_force_circle_fit = 0;
expected_radii_range = [1 10];
flag_use_permutations = [];

%% Find all the domains
fig_num = 234;
figure(fig_num); clf;
inputPoints = onearc_test_points;
domains_onearc_test_points  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (fig_num));

fig_num = 235;
figure(fig_num); clf;
inputPoints = corrupted_twoarc_test_points;
domains_corrupted_twoarc_test_points  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (fig_num));

% BASIC call with circle data, fitting it with a circle by not specifying station tolerance
fig_num = 1111;
figure(fig_num); clf;

inputPoints = corrupted_circle_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_use_permutations = [];


domains_corrupted_circle_test_points  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (fig_num));


%% Basic call with clean data
fig_num = 1;
figure(fig_num);
clf;
hold on;


regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_onearc_test_points{1}, fig_num);

true_params = [arc_true_circleCenter(1,1),arc_true_circleCenter(1,2), arc_true_circleRadius, arc_true_start_angle_in_radians, arc_true_end_angle_in_radians, 0.00 ];
fprintf(1,'\n\nResults of arc regression fitting:\n')
fprintf(1,'                        [centerX          centerY         radius         startAngle      endAngle        isCircle] (meters and degrees)\n');
params = true_params;
fcn_INTERNAL_printResults('True parameters', params)
params = domains_onearc_test_points{1}.best_fit_parameters;
fcn_INTERNAL_printResults('Hough parameters', params)
params = regression_domain.best_fit_parameters;
fcn_INTERNAL_printResults('Hough parameters', params)
fprintf(1,'ERRORS:\n');
params = abs( domains_onearc_test_points{1}.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Hough errors', params)
params = abs(regression_domain.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Regression errors', params)



%% Basic call with bad data
fig_num = 2;
figure(fig_num);
clf;
hold on;

regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_twoarc_test_points{1}, fig_num); 

true_params = [arc_true_circleCenter(1,1),arc_true_circleCenter(1,2), arc_true_circleRadius, arc_true_start_angle_in_radians, arc_true_end_angle_in_radians, 0.00 ];
fprintf(1,'\n\nResults of arc regression fitting:\n')
fprintf(1,'                        [centerX          centerY         radius         startAngle      endAngle        isCircle] (meters and degrees)\n');
params = true_params;
fcn_INTERNAL_printResults('True parameters', params)
params = domains_onearc_test_points{1}.best_fit_parameters;
fcn_INTERNAL_printResults('Hough parameters', params)
params = regression_domain.best_fit_parameters;
fcn_INTERNAL_printResults('Regression parameters', params)
fprintf(1,'ERRORS:\n');
params = abs( domains_onearc_test_points{1}.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Hough errors', params)
params = abs(regression_domain.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Regression errors', params)


%% Basic call with circle data
fig_num = 3;
figure(fig_num);
clf;
hold on;

regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_circle_test_points{1}, fig_num); 

true_params = circle_true_parameters(1,1:3);
fprintf(1,'\n\nResults of circle regression fitting:\n')
fprintf(1,'                        [centerX          centerY         radius         startAngle      endAngle        isCircle] (meters and degrees)\n');
params = true_params;
fcn_INTERNAL_printResults('True parameters', params)
params = domains_corrupted_circle_test_points{1}.best_fit_parameters;
fcn_INTERNAL_printResults('Hough parameters', params)
params = regression_domain.best_fit_parameters;
fcn_INTERNAL_printResults('Regression parameters', params)
fprintf(1,'ERRORS:\n');
params = abs( domains_corrupted_circle_test_points{1}.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Hough errors', params)
params = abs(regression_domain.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Regression errors', params)


%% Test of fast mode
% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;

    regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_twoarc_test_points{1}, fig_num); 

    
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;

    regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_twoarc_test_points{1}, fig_num); 

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitArcRegressionFromHoughFit:\n');
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

function fcn_INTERNAL_printResults(param_string, params)
header_string = sprintf('%s',param_string);
fixed_header_string = fcn_DebugTools_debugPrintStringToNCharacters(header_string,25);
fprintf(1,'%s ',fixed_header_string)
for ith_value = 1:length(params)
    param_string = sprintf('%.4f',params(ith_value));
    fixed_param_string = fcn_DebugTools_debugPrintStringToNCharacters(param_string,15);
    fprintf(1,'%s ',fixed_param_string)   
end
fprintf(1,'\n');
end