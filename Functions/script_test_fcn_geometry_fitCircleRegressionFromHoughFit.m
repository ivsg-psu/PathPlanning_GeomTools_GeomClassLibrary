% script_test_fcn_geometry_fitCircleRegressionFromHoughFit
% Exercises the function: fcn_geometry_fitCircleRegressionFromHoughFit

% Revision history:
% 2024_01_09 - S. Brennan
% -- wrote the code

close all;
clc;


%% Fill test data 
fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [3 4];
circle_radius = 2;
M = 50; % points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% Basic call with clean data
fig_num = 1;
figure(fig_num);
clf;
hold on;

[regression_fit_circle, domain_box, radial_errors, standard_deviation] = fcn_geometry_fitCircleRegressionFromHoughFit([circle_test_points(1,:); circle_test_points(2,:); circle_test_points(end,:)],circle_test_points, fig_num);

fprintf(1,'\n\nResults of circle regression fitting:\n')
fprintf(1,'Actual circle center [X Y] (meters): %.4f %.4f\n',circle_center(1,1),circle_center(1,2));
fprintf(1,'Fitted circle center [X Y] (meters): %.4f %.4f\n',regression_fit_circle(1,1),regression_fit_circle(1,2));
fprintf(1,'Predicted max distance error between actual and fitted center (meters)   %.4f\n',sum((circle_center - regression_fit_circle(1:2)).^2,2).^0.5);
fprintf(1,'Measured actual distance error between actual and fitted center (meters) %.4f\n',sigma/(length(circle_test_points(:,1))^0.5));
fprintf(1,'Actual circle radius (meters): %.4f \n',circle_radius);
fprintf(1,'Fitted circle radius (meters): %.4f \n',regression_fit_circle(1,3));
fprintf(1,'Radius distance error between actual and fitted (meters) %.4f\n',(circle_radius - regression_fit_circle(1,3)));

%% Basic call with bad data
fig_num = 2;
figure(fig_num);
clf;
hold on;

[regression_fit_circle, domain_box, radial_errors, standard_deviation] = fcn_geometry_fitCircleRegressionFromHoughFit([corrupted_circle_test_points(1,:); corrupted_circle_test_points(2,:); corrupted_circle_test_points(end,:)],corrupted_circle_test_points, fig_num);

fprintf(1,'\n\nResults of circle regression fitting:\n')
fprintf(1,'Actual circle center [X Y] (meters): %.4f %.4f\n',circle_center(1,1),circle_center(1,2));
fprintf(1,'Fitted circle center [X Y] (meters): %.4f %.4f\n',regression_fit_circle(1,1),regression_fit_circle(1,2));
fprintf(1,'Predicted max distance error between actual and fitted center (meters)   %.4f\n',sum((circle_center - regression_fit_circle(1:2)).^2,2).^0.5);
fprintf(1,'Measured actual distance error between actual and fitted center (meters) %.4f\n',sigma/(length(corrupted_circle_test_points(:,1))^0.5));
fprintf(1,'Actual circle radius (meters): %.4f \n',circle_radius);
fprintf(1,'Fitted circle radius (meters): %.4f \n',regression_fit_circle(1,3));
fprintf(1,'Radius distance error between actual and fitted (meters) %.4f\n',(circle_radius - regression_fit_circle(1,3)));



%% Test of fast mode
% Perform the calculation in slow mode
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [regression_fit_circle, domain_box, radial_errors, standard_deviation] = fcn_geometry_fitCircleRegressionFromHoughFit([circle_test_points(1,:); circle_test_points(2,:); circle_test_points(end,:)],circle_test_points, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [regression_fit_circle, domain_box, radial_errors, standard_deviation] = fcn_geometry_fitCircleRegressionFromHoughFit([circle_test_points(1,:); circle_test_points(2,:); circle_test_points(end,:)],circle_test_points, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitCircleRegressionFromHoughFit:\n');
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
