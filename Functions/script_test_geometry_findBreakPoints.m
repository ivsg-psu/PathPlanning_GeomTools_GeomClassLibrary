% script to find the break points


% Fill data points with lines and arcs
clc
% clear 
close all


% Assumptions: All the fits should have unique points. No input point
% should be repeated in any geometric fits.

rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

fig_num = 23;
 

% Arc 1 test points
% seed_points = [1 1; 2 1.8; 3 3];
seed_points = [1 1; 2.5 1.6; 3 3];
M = 5;
sigma = 0.02;

test_points_arc1 = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
hold on
% Line 1 test points
seed_points = [3 3; 6 4];
M = 10;
sigma = 0.02;

test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Arc 2 test points
seed_points = [6 4; 7.5 4.6; 8 6];
%seed_points = [2 3; 4 5; 6 3; 1 1];
M = 5;
sigma = 0.02;

test_points_arc2 = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_line1; test_points_arc1; test_points_arc2];

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);

kk = corrupted_testpoints(1:32,:);
ll = corrupted_testpoints(33:end,:);

corrupted_testpoints2 = [ll; kk];

%% Hough Segmentation
fig_num = 111;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


% q = domains{1,1}.points_in_domain;
% w = domains{1,2}.points_in_domain;
% e = domains{1,3}.points_in_domain;


%% 

% fucntion [breakPointsCell] = fcn_geometry_findBreakpoints(domains)

N_houghDomains = size(domains,2) - 1;

% Empty breakPoints domain structure
breakPoints.firstBreakPoint = [nan nan];
breakPoints.lastBreakPoint = [nan nan];
breakPoints.fitType = 'empty';
breakPoints.fitParameters = nan;

% Create a cell array to save the structure of each fit
breakPointsCell = cell(1,N_houghDomains);

% This loop saves the first and last points of each fit domain as the first
% and last break points in breakPointsCell cell array. This cell array also
% stores the fit type and parameters

for i = 1:N_houghDomains

    breakPointsCell{i} = breakPoints;
    breakPointsCell{i}.firstBreakPoint = domains{i}.points_in_domain(1,:);
    breakPointsCell{i}.lastBreakPoint = domains{i}.points_in_domain(end,:);
    breakPointsCell{i}.fitType = domains{i}.best_fit_type;
    breakPointsCell{i}.fitParameters = domains{i}.best_fit_parameters;

end


% end


