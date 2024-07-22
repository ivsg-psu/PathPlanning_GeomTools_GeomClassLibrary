% script_test_fcn_geometry_fitSegmentViaRegression
% Exercises the function: fcn_geometry_fitSegmentViaRegression

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

rng(1);

% Fill in points
seed_points = [0 0; 4 4];
M = 40;
sigma = 0.02;
line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);


% Test the regression
[best_fit_parameters, std_dev_orthogonal_distance] = fcn_geometry_fitSegmentViaRegression(line_test_points, fig_num);

% Check the output type and size
assert(isequal(size(best_fit_parameters),[1 4]));
assert(isequal(size(std_dev_orthogonal_distance),[1 1]));

assert(isequal(round(best_fit_parameters,1),round([0 0 pi/4 4*(2^0.5)],1)));
assert(isequal(round(std_dev_orthogonal_distance,1),round(sigma,1)));




%% Vertical line fit, down to up
fig_num = 3;
figure(fig_num); clf;

rng(1)

% Fill in points
seed_points = [2 3; 2 15];
M = 40;
sigma = 0.02;
line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

% Test the regression
[best_fit_parameters, std_dev_orthogonal_distance] = fcn_geometry_fitSegmentViaRegression(line_test_points, fig_num);

% Check the output type and size
assert(isequal(size(best_fit_parameters),[1 4]));
assert(isequal(size(std_dev_orthogonal_distance),[1 1]));

assert(isequal(round(best_fit_parameters,1),round([2 3 pi/2 12],1)));
assert(isequal(round(std_dev_orthogonal_distance,1),round(sigma,1)));


%% Vertical line fit, down to up
fig_num = 3;
figure(fig_num); clf;

rng(1)

% Fill in points
seed_points = [2 15; 2 3];
M = 40;
sigma = 0.02;
line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

% Test the regression
[best_fit_parameters, std_dev_orthogonal_distance] = fcn_geometry_fitSegmentViaRegression(line_test_points, fig_num);

% Check the output type and size
assert(isequal(size(best_fit_parameters),[1 4]));
assert(isequal(size(std_dev_orthogonal_distance),[1 1]));

assert(isequal(round(best_fit_parameters,1),round([2 15 3*pi/2 12],1)));
assert(isequal(round(std_dev_orthogonal_distance,1),round(sigma,1)));





%% Test of fast mode

rng(1)

% Fill in points
seed_points = [2 15; 2 3];
M = 40;
sigma = 0.02;
line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

% Perform the calculation in slow mode
fig_num = [];
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [best_fit_parameters, std_dev_orthogonal_distance] = fcn_geometry_fitSegmentViaRegression(line_test_points, fig_num);
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
     [best_fit_parameters, std_dev_orthogonal_distance] =  fcn_geometry_fitSegmentViaRegression(line_test_points, fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitSegmentViaRegression:\n');
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

