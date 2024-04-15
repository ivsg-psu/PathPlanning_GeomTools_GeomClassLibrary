% script_test_fcn_geometry_findEndPoints.m
% tests fcn_geometry_findEndPoints.m

% This script uses "fcn_geometry_HoughSegmentation.m" to fit the segments

% Revision History
% 2024_02_08 - Aneesh Batchu
% -- Started the script


% Note: All the fits should have unique points. No input point
% should be repeated in any geometric fits.
% In this case, the detected endpoints of the fits are valid. 


%% Clear workspace
close all

%% Two line segments: The ending point of a line segment is joined to the starting point of the other line segment 

rng(343)

fig_num = -1; 
% fig_num = 112;

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
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = -1;
% 
% figure(fig_num)


% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], fig_num);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

tolerance = 1; 

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, 113);

assert(isequal(round(closeEndPointsMatrix(2,:),4),[2.0933, 4.0370]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[1.9954, 3.9620]));
assert(isequal(round(dist_btw_close_endPoints,4), 0.1234));

%% Two line segments: The ending point of a line segment is not joined to the starting point of the other line segment 

rng(343)

fig_num = -1;
% fig_num = 115;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

tolerance = 1;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

% disp(endPointsCell)
% disp(closeEndPointsMatrix)
% disp(dist_btw_close_endPoints)

assert(isequal(round(closeEndPointsMatrix(2,:),4),[1.9954, 3.9620]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[2.4976, 4.5071]));
assert(isequal(round(dist_btw_close_endPoints,4), 0.7412));


%% Two line segments (Far): The ending point of a line segment is not joined to the starting point of the other line segment 

rng(343)

fig_num = -1;
% fig_num = 118;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


tolerance = 2.1;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

assert(isequal(round(closeEndPointsMatrix(1,:),4),[1.0351, 1.9824]));
assert(isequal(round(closeEndPointsMatrix(2,:),4),[6.0000, 5.0320]));
assert(isequal(round(dist_btw_close_endPoints,4), 2.2609));

%% Two line segments: Both the segments have an intersecting point

rng(343)

fig_num = -1;
% fig_num = 121;

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
fig_NuM = fig_num;

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

tolerance = 2;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

assert(isequal(round(closeEndPointsMatrix(1,:),4),[1.0042, 4.0062]));
assert(isequal(round(closeEndPointsMatrix(2,:),4),[1.9954, 3.9620]));
assert(isequal(round(dist_btw_close_endPoints,4), 2.8702));

%% Two line segments: Extend one of the line segments to find the intersecting point

rng(343)

fig_num = -1;
% fig_num = 124;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

tolerance = 1.1;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

assert(isequal(round(closeEndPointsMatrix(2,:),4),[1.9954, 3.9620]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[2.4958, 3.0062]));
assert(isequal(round(dist_btw_close_endPoints,4), 1.0789));

%% Two line segments: Both segments seem to be parallel but actually they are not (Slopes are not equal). Extend both the segments to find an intersecting point

rng(343)

fig_num = -1;
% fig_num = 127;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


tolerance = 2.5;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);


assert(isequal(round(closeEndPointsMatrix(1,:),4),[1.0373, 1.9876]));
assert(isequal(round(closeEndPointsMatrix(2,:),4),[4.0045, 5.9895]));
assert(isequal(round(dist_btw_close_endPoints,4), 2.5089));

%% Two line segments: Both segments are parallel (Slopes are equal).

rng(343)

fig_num = -1;
% fig_num = 130;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

tolerance = 3.2;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

assert(isequal(round(closeEndPointsMatrix(2,:),4),[4.0000, 1.9768]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[5.0000, 4.9955]));
assert(isequal(round(dist_btw_close_endPoints,4), 3.1800));

%% Two line segments: Both segments are parallel (Slopes are equal) and close to each other

rng(343)

fig_num = -1;
% fig_num = 133;

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
fig_NuM = fig_num;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),fig_num);


% Hough Segmentation
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


tolerance = 0.5;
% [endPointsCell, closeEndPointPairs] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

[~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

assert(isequal(round(closeEndPointsMatrix(2,:),4),[4.0000, 1.9768]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[4.0000, 2.1955]));
assert(isequal(round(dist_btw_close_endPoints,4), 0.2187));

%% (Fail CASE) Four Line segements (Double yellow): Two pairs of parallel segments are close to each other but not intersecting

if 1 == 0
  
    fig_num = -1;
    % fig_num = 136;

    % Line 1 test points - Pair 1
    seed_points = [1 2; 4 2];
    M = 10;
    sigma = 0.02;

    test_points_LineSegment1_pair1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


    % Line 2 test points - Pair 1
    seed_points = [1 2.2; 4 2.2];
    M = 10;
    sigma = 0.02;

    test_points_LineSegment2_pair1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


    % Line 3 test points - Pair 2
    seed_points = [5 2; 8 2];
    M = 10;
    sigma = 0.02;

    test_points_LineSegment3_pair2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


    % Line 4 test points - Pair 2
    seed_points = [5 2.2; 8 2.2];
    M = 10;
    sigma = 0.02;

    test_points_LineSegment4_pair2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

    testpoints = [test_points_LineSegment1_pair1; test_points_LineSegment2_pair1; test_points_LineSegment3_pair2; test_points_LineSegment4_pair2];

    % Corrupt the points
    fig_num = 137;
    fig_NuM = fig_num;
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 3;

    corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
        (probability_of_corruption), (magnitude_of_corruption),fig_num);

    % Hough Segmentation
    fig_num = -1;
    transverse_tolerance = 0.05; % Units are meters
    station_tolerance = 0.5; % Units are meters.
    threshold_max_points = 10;
    input_points = corrupted_testpoints;

    domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);



    tolerance = 1.2;
    % [endPointsCell, closeEndPointPairs, distance_btw_breakpoints] = fcn_geometry_findEndPoints(domains, tolerance, fig_num);

    [~, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, [], -1);

    endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

    endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
    endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

    [closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, fig_NuM);

    % assert(isequal(round(closeEndPointsMatrix(2,:),4),[4.0000, 1.9768]));
    % assert(isequal(round(closeEndPointsMatrix(3,:),4),[4.0000, 2.1955]));
    % assert(isequal(round(dist_btw_close_endPoints,4), 0.2187));

end
%% In this case, one arc and one line have been used as the test data. 

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

fig_num = -1;


% Arc 1 test points
% seed_points = [1 1; 2 1.8; 3 3];
seed_points = [2 3; 4 5; 6 3];
M = 10; % Points per meter
sigma = 0.02;

% Fill test data for arc 1
[test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

fig_num = 1112;
corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Line 1 test points
seed_points = [5 6; 9 3];
M = 10;
sigma = 0.02;

fig_num = -1;
% Fill test data for line 1
test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


fig_num = 1112;
corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);


% Hough Segmentation
fig_num = -1;
% fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1.5; % Units are meters. 
threshold_max_points = 20;
input_points = [corrupted_test_points_arc1; corrupted_test_points_line1];

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = -1;
% fig_num = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], fig_num);
fcn_geometry_plotFitDomains(regression_domains, -1);

fig_num = -1;
tolerance = 3.2;

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(regression_domains, [], fig_num);

endPointsMatrix = zeros(numel(sortedHoughSegmentEndPoints)/2, 2);

endPointsMatrix(1:2:end,:) = sortedHoughSegmentEndPoints(:,1:2);
endPointsMatrix(2:2:end,:) = sortedHoughSegmentEndPoints(:,3:4);

[closeEndPointsMatrix, dist_btw_close_endPoints] = fcn_geometry_findEndPoints(endPointsMatrix(1,:), endPointsMatrix(2:end-1,:), endPointsMatrix(end,:), tolerance, 1112);


assert(isequal(round(closeEndPointsMatrix(2,:),4),[6.0106, 3.0000]));
assert(isequal(round(closeEndPointsMatrix(3,:),4),[ 5.0026, 6.0035]));
assert(isequal(round(dist_btw_close_endPoints,4), 3.1681));



% 
% %%
% function print_results(distance,location)
% fprintf(1,'Distance \t Location X \t Location Y \n');
% if ~isempty(distance)
%     for i_result = 1:length(distance(:,1))
%         fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
%     end
% end
% end
% 
% %%
% function print_more_results(distance,location,path_segments)
% fprintf(1,'Distance \t Location X \t Location Y \t PathSegment \n');
% if ~isempty(distance)
%     for i_result = 1:length(distance(:,1))
%         fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
%     end
% end
% end