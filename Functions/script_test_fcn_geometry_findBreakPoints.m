% script_test_fcn_geometry_findBreakPoints.m
% tests fcn_geometry_findBreakPoints.m

% Revision History
% 2024_02_08 - Aneesh Batchu
% -- Started the script


% Note: All the fits should have unique points. No input point
% should be repeated in any geometric fits.
% In this case, the detected breakpoints of the fits are valid. 


%% Clear workspace
clc
clear 
close all

%% Simple Line Case

rng(343)

fig_num = 13;

% Line test points
seed_points = [1 2; 4 4];
M = 10;
sigma = 0.02;

test_points_line = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the points
fig_num = 14;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(test_points_line,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 111;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

% Find break points

breakPointsCell = fcn_geometry_findBreakpoints(domains);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

%% Simple Arc Case

rng(343)

fig_num = 23;

seed_points = [1 1; 2.5 1.6; 3 3];
M = 5;
sigma = 0.02;

test_points_arc = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the points
fig_num = 24;
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(test_points_arc,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 222;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

% Find break points

breakPointsCell = fcn_geometry_findBreakpoints(domains);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

%% In this case, two arcs and one line have been used as the test data. The break points are saved in breakPointsCell array 

% Fill data points with lines and arcs
rng(3)

probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

fig_num = 33;
 

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

% corrupt the test points
fig_num = 34;
corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);

% kk = corrupted_testpoints(1:32,:);
% ll = corrupted_testpoints(33:end,:);
% 
% corrupted_testpoints2 = [ll; kk];

% Hough Segmentation
fig_num = 333;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


% q = domains{1,1}.points_in_domain;
% w = domains{1,2}.points_in_domain;
% e = domains{1,3}.points_in_domain;

[breakPointsCell, breakPointClosePairs] = fcn_geometry_findBreakpoints(domains);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(domains{2}.points_in_domain(1,:) == breakPointsCell{2}.firstBreakPoint, [1 1]));
assert(isequal(domains{2}.points_in_domain(end,:) == breakPointsCell{2}.lastBreakPoint, [1 1]));

assert(isequal(domains{3}.points_in_domain(1,:) == breakPointsCell{3}.firstBreakPoint, [1 1]));
assert(isequal(domains{3}.points_in_domain(end,:) == breakPointsCell{3}.lastBreakPoint, [1 1]));


%% More complicated case: Two arcs and Two straight lines (not perfect) 


% Hough Segmentation is not executed correctly. Need to give
% station-transverse coordinates as the input points to get the unique fit
% in each stage. 
% However, the findBreakPoints function finds the break points. The break
% points are not valid since the geometric fits are not unique.

rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

fig_num = 43;
 

% Arc 1 test points
% seed_points = [1 1; 2 1.8; 3 3];
seed_points = [3 1; 1.5 2.5; 3 4];
M = 10;
sigma = 0.02;

test_points_arc1 = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

hold on

% Line 1 test points
seed_points = [3 4; 7 3];
M = 10;
sigma = 0.02;

test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Line 2 test points
seed_points = [3 1; 7 2];
M = 10;
sigma = 0.02;

test_points_line2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Arc 2 test points
seed_points = [7 2; 7.5 2.5; 7 3];

M = 10;
sigma = 0.02;

test_points_arc2 = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_line1; test_points_line2; test_points_arc1; test_points_arc2];

fig_num = 44;
corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);

% Hough Segmentation
fig_num = 444;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1; % Units are meters. 
threshold_max_points = 5;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

breakPointsCell = fcn_geometry_findBreakpoints(domains);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(domains{2}.points_in_domain(1,:) == breakPointsCell{2}.firstBreakPoint, [1 1]));
assert(isequal(domains{2}.points_in_domain(end,:) == breakPointsCell{2}.lastBreakPoint, [1 1]));

assert(isequal(domains{3}.points_in_domain(1,:) == breakPointsCell{3}.firstBreakPoint, [1 1]));
assert(isequal(domains{3}.points_in_domain(end,:) == breakPointsCell{3}.lastBreakPoint, [1 1]));

assert(isequal(domains{4}.points_in_domain(1,:) == breakPointsCell{4}.firstBreakPoint, [1 1]));
assert(isequal(domains{4}.points_in_domain(end,:) == breakPointsCell{4}.lastBreakPoint, [1 1]));
