% script_test_fcn_geometry_fitLinearRegressionFromHoughFit
% Exercises the function: fcn_geometry_fitLinearRegressionFromHoughFit

% Revision history:
% 2023_12_15 - S. Brennan
% -- wrote the code
% 2024_04_11 - S. Brennan
% -- added assertion testing

close all;


%% Basic call - line fitting

fig_num = 1;
figure(fig_num);
clf;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Use Hough line fitting to separate inliners and outliers
transverse_tolerance = 0.05;
station_tolerance = 2;
points_required_for_agreement = 20;

domains_line_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);


regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_line_fitting{1}, [], fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_line_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));


%% Show no figure is generated

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Use Hough line fitting to separate inliners and outliers
transverse_tolerance = 0.05;
station_tolerance = 2;
points_required_for_agreement = 20;

domains_line_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);


regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_line_fitting{1});

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_line_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));

%% Now plot the results

fig_num = 2;
figure(fig_num);
clf;
hold on;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Use Hough line fitting to separate inliners and outliers

transverse_tolerance = 0.05;
station_tolerance = 2;
points_required_for_agreement = 20;

domains_line_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);


regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_line_fitting{1}, [], fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_line_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));


%% Vertical line fit

fig_num = 3;
figure(fig_num); clf;

% Fill in points
seed_points = [2 3; 2 15];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Use Hough line fitting to separate inliners and outliers

transverse_tolerance = 0.05;
station_tolerance = 2;
points_required_for_agreement = 20;

domains_line_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);


[regression_domain, std_dev_transverse_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(domains_line_fitting{1}, [], fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_line_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));

fprintf(1,'\n\nFitting results: \n');
fprintf(1,'Expected standard deviation in fit, transverse direction (total least squares), in meters: %.4f\n',sigma);
fprintf(1,'Measured standard deviation in fit, transverse direction (total least squares), in meters: %.4f\n',std_dev_transverse_distance);

%% Basic call - segment fitting

% Fill test data 
fig_num = 4;
figure(fig_num);
clf;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Demo Hough line segment fitting

transverse_tolerance = 0.1;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains_segment_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);

regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_segment_fitting{1}, [], fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_segment_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));

%% Test the best_fit_domain_box_projection_distance input

% Fill test data 
fig_num = 5;
figure(fig_num);
clf;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Demo Hough line segment fitting

transverse_tolerance = 0.1;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains_segment_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);

best_fit_domain_box_projection_distance = 2; % A HUGE number
regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_segment_fitting{1}, best_fit_domain_box_projection_distance, fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_segment_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));

%% Fill test data 

fig_num = 6;
figure(fig_num);
clf;
hold on;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Demo Hough line segment fitting

transverse_tolerance = 0.1;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains_segment_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);

regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_segment_fitting{1}, [], fig_num);

% Check the output type and size
assert(isstruct(regression_domain));
assert(isfield(regression_domain,'best_fit_type'));
assert(isfield(regression_domain,'points_in_domain'));
assert(isfield(regression_domain,'best_fit_parameters'));
assert(isfield(regression_domain,'best_fit_domain_box'));
assert(isfield(regression_domain,'best_fit_1_sigma_box'));
assert(isfield(regression_domain,'best_fit_2_sigma_box'));
assert(isfield(regression_domain,'best_fit_3_sigma_box'));
assert(ischar(regression_domain.best_fit_type));
assert(isequal(regression_domain.points_in_domain,domains_segment_fitting{1}.points_in_domain));
assert(isequal(size(regression_domain.best_fit_parameters),[1 4]));
assert(issimplified(regression_domain.best_fit_domain_box));
assert(issimplified(regression_domain.best_fit_1_sigma_box));
assert(issimplified(regression_domain.best_fit_2_sigma_box));
assert(issimplified(regression_domain.best_fit_3_sigma_box));


%% Test of fast mode

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Demo Hough line segment fitting

transverse_tolerance = 0.1;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains_segment_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, -1);


% Perform the calculation in slow mode
fig_num = [];
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [regression_domain, std_dev_transverse_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(domains_segment_fitting{3}, [], fig_num);
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
    [regression_domain, std_dev_transverse_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(domains_segment_fitting{3}, [], fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitLinearRegressionFromHoughFit:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

% %% Debugging example
% % Noisy vertical line
% fig_num = 999;
% points = [
%    9.991259411578580   0.200000000000000
%   10.001803670408686   0.300000000000000
%    9.982375541376497   0.400000000000000
%    9.977758789496919   0.500000000000000
%   10.008682631868705   0.700000000000000
%    9.988763772969360   0.800000000000000
%   10.003402293646682   1.000000000000000
%   10.011493908019808   1.200000000000000
%    9.984955058249714   1.300000000000000
%    9.995849995905278   1.400000000000000
%   10.006337144959792   1.500000000000000
%    9.982672092520335   1.600000000000000
%   10.007135029591263   1.700000000000000
%    9.998450373432329   1.800000000000000
%   10.005582924433416   1.900000000000000
%    9.979543989497849   2.000000000000000
%    9.974383738360030   2.100000000000000
%   10.001778511906023   2.500000000000000
%    9.988736840713367   2.600000000000000
%    9.989763064150692   2.800000000000000
%   10.009874363357989   3.000000000000000
%    9.984802416100781   3.100000000000000
%    9.998229855546459   3.300000000000000
%    9.984344000425939   3.400000000000000
%    9.975609845675763   3.500000000000000
%    9.991928970032030   3.600000000000000
%    9.990306079204805   3.900000000000000
%    9.983588133844744   4.000000000000000
%   10.005622755909043   4.100000000000000
%    9.981054848371585   4.200000000000000
%    9.993381147363705   4.300000000000000
%   10.007368013343468   4.400000000000000
%   10.007946121327633   4.500000000000000
%    9.990592139085731   4.600000000000000
%   10.009923513002883   4.700000000000000
%    9.971033609842408   4.900000000000000
%     ];
% 
% 
% % Create dummy data
% test_domain = fcn_geometry_fillEmptyDomainStructure;
% test_domain.best_fit_type = 'Hough line';
% test_domain.points_in_domain = points;
% test_domain.best_fit_parameters = [10 0 10 5];
% 
% 
% [regression_domain, std_dev_transverse_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(test_domain, fig_num);
% fprintf(1,'\n\nFitting results: \n');
% fprintf(1,'Standard deviation in fit, transverse direction (total least squares), in meters: %.4f\n',std_dev_transverse_distance);
% 


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

