% script_test_fcn_geometry_calcLinearRegressionFromHoughFit
% Exercises the function: fcn_geometry_calcLinearRegressionFromHoughFit

% Revision history:
% 2023_12_15
% -- wrote the code

close all;
clc;


%% Fill test data - 1 segment
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5];
M = 20; % M is points per meter
sigma = 0.1;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,fig_num);

% Add outliers to corrupt the results
probability_of_corruption = 0.05;
magnitude_of_corruption = 1;

test_points_with_outliers = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

%% Basic call 
fig_num = 1;
figure(fig_num);
clf;
hold on;

plot(test_points_with_outliers(:,1),test_points_with_outliers(:,2),'k.','MarkerSize',20);
[regression_fit_line_segment, domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit([test_points(1,:); test_points(end,:)],test_points, fig_num);

%% Test of fast mode
% Perform the calculation in slow mode
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [regression_fit_line_segment, domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit([test_points(1,:); test_points(end,:)],test_points, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [regression_fit_line_segment, domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit([test_points(1,:); test_points(end,:)],test_points, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_calcLinearRegressionFromHoughFit:\n');
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

