% script_test_fcn_geometry_plotFitDomains
% Exercises the function: fcn_geometry_plotFitDomains

% Revision history:
% 2024_01_12
% -- wrote the code

close all;
clc;

%% Fill in some test data
% This takes a while - it's generated from the test script for Hough
% Segmentation

% clear example_domains

if ~exist('example_domains','var')
    % Advanced example 3: find segments within a chevron
    M = 10; % 40 points per meter

    rng(234)
    sigma = 0.02;

    multi_segment_test_points = [];

    seed_points = [0 0; 10 0];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points]; %#ok<*NASGU>

    seed_points = [0 0; 10 5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [2 0; 3 1.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [4 0; 5 2.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [6 0; 7 3.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [8 0; 9 4.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [10 0; 10 5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];
     
    % Add outliers to corrupt the results
    outliers = [10*randn(100,1) 5*randn(100,1)];
    multi_segment_test_points = [multi_segment_test_points; outliers];


    % Call the segmentation function
    fig_num = 3;
    transverse_tolerance = 0.1; % Units are meters
    station_tolerance = 0.2; % Units are meters
    threshold_max_points = 10;
    input_points = multi_segment_test_points;

    example_domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);
end

%% BASIC test of plotting
fig_num = 1234;
figure(fig_num);
clf;
hold on;
grid on;
axis equal
grid minor;

plot(multi_segment_test_points(:,1),multi_segment_test_points(:,2),'k.','MarkerSize',20);
fcn_geometry_plotFitDomains(example_domains, fig_num);




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

