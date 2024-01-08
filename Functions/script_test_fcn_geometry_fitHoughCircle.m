% script_test_fcn_geometry_circleHoughFit
% Exercises the function: fcn_geometry_circleHoughFit

% Revision history:
% 2023_12_15
% -- wrote the code
% 2024_01_04
% -- fixed the argument inputs

%% Set up the workspace

clc
close all

%% Examples for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ______                           _
% |  ____|                         | |
% | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                            | |
%                            |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Fill test data 
fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 5; % 5 points per meter
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

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Fill test data - 2 arcs
corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];

%% Example - 1 - BASIC call with arc data, fitting it with a circle
fig_num = 111;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
expected_radii_range = [];
flag_use_permutations = [];

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, flag_is_a_circle] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

%% Example - 1 - BASIC call with circle data, fitting it with a circle
fig_num = 222;
figure(fig_num); clf;

inputPoints = corrupted_circle_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
expected_radii_range = [];
flag_use_permutations = [];

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);


%% Test using expected radii range
fig_num = 2;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
expected_radii_range = [1 3];
flag_use_permutations = [];

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

%% Speed test effect of adding radii range to show this speeds up calculations
% Note, there is more speed-up the more corrupted and larger the data is
% used
inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;
expected_radii_range = [1 3];
flag_use_permutations = [];

% Perform the calculation in slow mode
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, [], flag_use_permutations, -1);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fitHoughCircle without radii range (slow) and with radii range (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Test arc fitting using station_tolerance

inputPoints = corrupted_twoarc_test_points;
transverse_tolerance = 0.1;
expected_radii_range = [1 3];
flag_use_permutations = [];

% Use station tolerance low to find only largest arc
station_tolerance = 0.3;
fig_num = 7777;
figure(fig_num); clf;
[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

% Make station tolerance larger so it finds entire arc, connecting together
% but not finding a circle
station_tolerance = 3;
fig_num = 7788;
figure(fig_num); clf;
[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

% Fit a circle by shutting station tolerance off
station_tolerance = [];
fig_num = 7799;
figure(fig_num); clf;
[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);



%% Test using flag_use_permutations to search only points in sequence
% This forces the search to assume the points are ordered. This speeds up
% the search process quite a bit, since much of the search is simply
% sorting points. This will not work well if the points are not well
% sorted, for example if there are a large number of outliers.

fig_num = 44;
figure(fig_num); clf;

inputPoints = onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 0.5;
expected_radii_range = [1 3];
flag_use_permutations = 0;

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

%% Test using flag_use_permutations for fractional setting (e.g. only 50%)
% This assumes that the points are over-fitted, e.g. that there are way
% more points than needed to fit. By setting a fraction, we can specify to
% only use a fraction of the input data, which speeds things up.

fig_num = 444;
figure(fig_num); clf;

inputPoints = onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 1;
expected_radii_range = [1 3];
flag_use_permutations = 0.5;

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);
fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters(1),trueParameters(2),trueParameters(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',best_fitted_parameters(1),best_fitted_parameters(2),best_fitted_parameters(3));

%% Test using flag_use_permutations for N-point setting
% This specifies the number of points to use as N points. This should only
% be used if the data is VERY clean, wherein all the data is quite good.
% However, if this is the case, this is quite fast.

fig_num = 55;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
expected_radii_range = [1 3];
flag_use_permutations = 30;

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);
fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters(1),trueParameters(2),trueParameters(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',best_fitted_parameters(1),best_fitted_parameters(2),best_fitted_parameters(3));

% Now show how it works with "clean" points
fig_num = 66;
figure(fig_num); clf;

inputPoints = onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
expected_radii_range = [1 3];
flag_use_permutations = 30;

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);
fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters(1),trueParameters(2),trueParameters(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',best_fitted_parameters(1),best_fitted_parameters(2),best_fitted_parameters(3));





%% Speed test effect of flag_use_permutations
seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 5; % Number of points per meter
sigma = 0.02;

corrupted_onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(corrupted_onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;
expected_radii_range = [1 3];
slow_flag_use_permutations = 1;
fast_flag_use_permutations = 0;

% Perform the calculation in slow mode
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, slow_flag_use_permutations, -1);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
slowParameters = best_fitted_parameters;
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 3; minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, fast_flag_use_permutations, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
fastParameters = best_fitted_parameters;
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fitHoughCircle with flag_use_permutations = %.2d (slow) and with flag_use_permutations=%.2f (fast):\n',slow_flag_use_permutations,fast_flag_use_permutations);
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'True parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters(1),trueParameters(2),trueParameters(3));
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Slow parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',slowParameters(1),slowParameters(2),slowParameters(3));
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Fast parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',fastParameters(1),fastParameters(2),fastParameters(3));
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Test of fast implementation mode 
inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;

% Perform the calculation in slow mode
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range,  flag_use_permutations, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 3; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitHoughCircle:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

