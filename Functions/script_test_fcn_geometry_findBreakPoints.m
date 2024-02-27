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
close all

%% Simple Line Segment Case

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
tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(closeBreakPointPairs,  zeros(0, 2)))


%% Two line segments: The ending point of a line segment is joined to the starting point of the other line segment 

rng(343)

fig_num = 112;

% Line test points
seed_points = [1 2; 2 4; 4 5];
M = 10;
sigma = 0.02;

test_points_twoLineSegments = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Corrupt the points
fig_num = 113;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(test_points_twoLineSegments,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 114;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 113;

figure(fig_num)

tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);

figure(fig_num)
hold on
plot(breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)


fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.firstBreakPoint(1) breakPointsCell{1}.firstBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 233;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


% Fit parameters of Line Segment
%
% breakPointsCell{}.fitParameters(1) = unit_projection_vector_x 
% breakPointsCell{}.fitParameters(2) = unit_projection_vector_y 
% breakPointsCell{}.fitParameters(3) = base_point_x 
% breakPointsCell{}.fitParameters(4) = base_point_y
% breakPointsCell{}.fitParameters(5) = station_distance_min 
% breakPointsCell{}.fitParameters(6) = station_distance_max

% Creating a matrix to store the best fit parameters
fitParametersMatrix = zeros([length(breakPointsCell),6]);

for i = 1: length(breakPointsCell)
    fitParametersMatrix(i,:) = breakPointsCell{i}.fitParameters;
    
end

% Finding the slopes of the lines
slopesLines = fitParametersMatrix(:,2)./fitParametersMatrix(:,1);
% Finding the y-intercepts of the lines 
yInterceptsLines = fitParametersMatrix(:,4) - slopesLines.*fitParametersMatrix(:,3);

% Line: ax+by+c = 0 ---> -mx+y-yIntercept = 0
% 
% a = -slopesLines; b = 1; c = -yInterceptsLines
%
% -mx + b = c
% AA = [-m1  1; -m2 1]
% bb = [c1; c2]

% Solving for intersecting point 
AA = [-1.*(slopesLines), ones(size(slopesLines))]; 
bb = (yInterceptsLines);

% Finding the intersection point
intersection_point = AA\bb;

% Printing the point
fprintf('Intersection Point (using backslash operator): (%.4f, %.4f)\n', intersection_point(1), intersection_point(2));

% Plotting the intersection point using the best fit parameters
figure(fig_num)
hold on 

% Define x range for plotting
x_range = linspace(-10, 10, 100);

% Calculate y values for each line
y1 = slopesLines(1,1) * x_range + yInterceptsLines(1,1);
y2 = slopesLines(2,1) * x_range + yInterceptsLines(2,1);

% Plot the lines
% figure(888)
plot(x_range, y1, 'k--', 'LineWidth', 2);
plot(x_range, y2, 'b--', 'LineWidth', 2); % plot second line in red

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)

%% Two line segments: The ending point of a line segment is not joined to the starting point of the other line segment 

rng(343)

fig_num = 115;

% Line 1 test points
seed_points = [1 2; 2 4];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [2.5 4.5; 4 5];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 116;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 117;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 116;
tolerance = 1;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);

figure(fig_num)
hold on
plot(breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)


intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

% Plotting the intersection point using the best fit parameters
figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)



fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.firstBreakPoint(1) breakPointsCell{1}.firstBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 333;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Two line segments (Far): The ending point of a line segment is not joined to the starting point of the other line segment 

rng(343)

fig_num = 118;

% Line 1 test points
seed_points = [1 2; 2 4];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [4 5; 6 5];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 119;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 120;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 119;
tolerance = 2.1;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);


figure(fig_num)
hold on
plot(breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)



intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)


fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.firstBreakPoint(1) breakPointsCell{1}.firstBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 433;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Two line segments: Both the segments have an intersecting point

rng(343)

fig_num = 121;

% Line 1 test points
seed_points = [1 2; 2 4];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [1 4; 4 2];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 122;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 123;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 122;

tolerance = 1;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);

figure(fig_num)
hold on
plot(breakPointsCell{2}.firstBreakPoint(1), breakPointsCell{2}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{1}.lastBreakPoint(1), breakPointsCell{1}.lastBreakPoint(2), 'b.', 'MarkerSize',20)


intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 
% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)




fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.lastBreakPoint(1), breakPointsCell{1}.lastBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.firstBreakPoint(1), breakPointsCell{2}.firstBreakPoint(2)];
fig_debugging = 533;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Two line segments: Extend one of the line segments to find the intersecting point

rng(343)

fig_num = 124;

% Line 1 test points
seed_points = [1 2; 2 4];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [2.5 3; 4 4];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 125;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 126;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 125;
tolerance = 1.1;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);


figure(fig_num)
hold on
plot(breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)



intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)




fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 633;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Two line segments: Both segments seem to be parallel but actually they are not (Slopes are not equal). Extend both the segments to find an intersecting point

rng(343)

fig_num = 127;

% Line 1 test points
seed_points = [1 2; 2 5];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [2.5 2.5; 4 6];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 128;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 129;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 128;
tolerance = 1.6;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);

figure(fig_num)
hold on
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{1}.lastBreakPoint(1), breakPointsCell{1}.lastBreakPoint(2), 'b.', 'MarkerSize',20)


intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)



fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.lastBreakPoint(1), breakPointsCell{1}.lastBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 733;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Two line segments: Both segments are parallel (Slopes are equal).

rng(343)

fig_num = 130;

% Line 1 test points
seed_points = [1 2; 4 2];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [5 5; 8 5];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 131;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 132;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 131;
tolerance = 3.1;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);


intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)

% figure(fig_num)
% hold on
% plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'c.', 'MarkerSize',20)
% plot(breakPointsCell{1}.lastBreakPoint(1), breakPointsCell{1}.lastBreakPoint(2), 'c.', 'MarkerSize',20)

fprintf(1,'Intersection result: \n');
wall_start = closeBreakPointPairs(1,:);
wall_end   = closeBreakPointPairs(2,:);
sensor_vector_start = closeBreakPointPairs(3,:);
sensor_vector_end   = closeBreakPointPairs(4,:);
fig_debugging = 833;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Two line segments: Both segments are parallel (Slopes are equal) and close to each other

rng(343)

fig_num = 133;

% Line 1 test points
seed_points = [1 2; 4 2];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [4 2.2; 7 2.2];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
fig_num = 134;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 135;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 134;
tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance, fig_num);

figure(fig_num)
hold on
plot(breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2), 'b.', 'MarkerSize',20)
plot(breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2), 'b.', 'MarkerSize',20)


intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)


fprintf(1,'Intersection result: \n');
wall_start = [breakPointsCell{1}.firstBreakPoint(1), breakPointsCell{1}.firstBreakPoint(2)];
wall_end   = closeBreakPointPairs(1,:);
sensor_vector_start = closeBreakPointPairs(2,:);
sensor_vector_end   = [breakPointsCell{2}.lastBreakPoint(1), breakPointsCell{2}.lastBreakPoint(2)];
fig_debugging = 933;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Three line segments: All three segments are close to each other

rng(343)

fig_num = 133;

% Line 1 test points
seed_points = [1 2; 4 2];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Line 2 test points
seed_points = [3 2.5; 6 2.5];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

% Line 3 test points
seed_points = [2 3; 8 3];
M = 10;
sigma = 0.02;

test_points_LineSegment3 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

testpoints = [test_points_LineSegment1; test_points_LineSegment2; test_points_LineSegment3];

% Corrupt the points
fig_num = 134;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = 135;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

tolerance = 0.5;
[breakPointsCell, ~] = fcn_geometry_findBreakpoints(domains, tolerance);

fig_num = 134;
intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num); 

figure(fig_num)
hold on 

% Plot the intersection point
plot(intersection_point(1), intersection_point(2), 'go', 'MarkerSize',30, 'LineWidth',2)
plot(intersection_point(1), intersection_point(2), 'c.', 'MarkerSize',20)

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

tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(closeBreakPointPairs,  zeros(0, 2)))

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

tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(domains{2}.points_in_domain(1,:) == breakPointsCell{2}.firstBreakPoint, [1 1]));
assert(isequal(domains{2}.points_in_domain(end,:) == breakPointsCell{2}.lastBreakPoint, [1 1]));

assert(isequal(domains{3}.points_in_domain(1,:) == breakPointsCell{3}.firstBreakPoint, [1 1]));
assert(isequal(domains{3}.points_in_domain(end,:) == breakPointsCell{3}.lastBreakPoint, [1 1]));


figure(34)

plot(closeBreakPointPairs(:,1), closeBreakPointPairs(:,2), '.c','MarkerSize',13);


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

tolerance = 0.5;
[breakPointsCell, closeBreakPointPairs] = fcn_geometry_findBreakpoints(domains, tolerance);

assert(isequal(domains{1}.points_in_domain(1,:) == breakPointsCell{1}.firstBreakPoint, [1 1]));
assert(isequal(domains{1}.points_in_domain(end,:) == breakPointsCell{1}.lastBreakPoint, [1 1]));

assert(isequal(domains{2}.points_in_domain(1,:) == breakPointsCell{2}.firstBreakPoint, [1 1]));
assert(isequal(domains{2}.points_in_domain(end,:) == breakPointsCell{2}.lastBreakPoint, [1 1]));

assert(isequal(domains{3}.points_in_domain(1,:) == breakPointsCell{3}.firstBreakPoint, [1 1]));
assert(isequal(domains{3}.points_in_domain(end,:) == breakPointsCell{3}.lastBreakPoint, [1 1]));

assert(isequal(domains{4}.points_in_domain(1,:) == breakPointsCell{4}.firstBreakPoint, [1 1]));
assert(isequal(domains{4}.points_in_domain(end,:) == breakPointsCell{4}.lastBreakPoint, [1 1]));

figure(44)
plot(closeBreakPointPairs(:,1), closeBreakPointPairs(:,2), '.c','MarkerSize',13);


%%
function print_results(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end

%%
function print_more_results(distance,location,path_segments)
fprintf(1,'Distance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end