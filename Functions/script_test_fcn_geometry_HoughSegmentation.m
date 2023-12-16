% script_test_fcn_geometry_HoughSegmentation
% Exercises the function: fcn_geometry_HoughSegmentation

% Revision history:
% 2023_12_15
% -- wrote the code

close all;
clc;


%% Fill test data - 3 segments
% fig_num = 23;
% figure(fig_num);
% clf;
% hold on;
% axis equal
% grid on;

seed_points = [2 3; 4 5; 7 0; 9 5; 10 20; 13 14];
M = 10; % 10 points per meter
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

%% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));

%% Basic call 
fig_num = 1;
transverse_tolerance = 0.02; % Units are meters
station_tolerance = 0.5; % Units are meters
threshold_max_points = 10;
input_points = test_points;
fig_points_plot = 11111;
domains = fcn_geometry_HoughSegmentation(test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


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

