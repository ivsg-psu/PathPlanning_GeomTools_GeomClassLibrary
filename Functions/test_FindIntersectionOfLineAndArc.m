% 
% 
% This script is written to develop a function that can find the
% intersection of an arc with a segment. 

% Revision History
% 2024_02_08 - Aneesh Batchu
% -- Started the script

%% clear workspace
clc
close all

%% In this case, one arc and one line have been used as the test data. 

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% fig_num = 111;
fig_num = -1;

% Arc 1 test points
% seed_points = [1 1; 2 1.8; 3 3];
seed_points = [2 3; 4 5; 6 3];
M = 10; % Points per meter
sigma = 0.02;

% Fill test data for arc 1
[test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

fig_num = 112;
corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Line 1 test points
seed_points = [5 6; 9 3];
M = 10;
sigma = 0.02;

% fig_num = 111;
fig_num = -1;
% Fill test data for line 1
test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


fig_num = 112;
corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);


% Hough Segmentation
% fig_num = 501;
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1.5; % Units are meters. 
threshold_max_points = 20;
input_points = [corrupted_test_points_arc1; corrupted_test_points_line1];

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 601;
% fig_num = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, fig_num);
% fcn_geometry_plotFitDomains(regression_domains, fig_num+2);

fig_num = 112;
tolerance = [];

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortHoughSegments(regression_domains, tolerance, fig_num);


hold on 

% Location of S (Start Point of Line Segment) 
pointS = endPointsCell{2}.firstEndPoint;

% Location of E (End Point of Line Segement)
pointE = endPointsCell{2}.lastEndPoint;

% Location of Arc center
pointC = endPointsCell{1}.fitParameters(1,1:2);

% plot Arc Center
plot(pointC(1,1), pointC(1,2), 'r.','MarkerSize',30)

% Calculate vector EC
vectorEC = pointC - pointE;

% Calculate vector ES
vectorES = pointS - pointE;

% vector_magnitude = sum(input_vector.^2,2).^0.5;
% unitVector_input_vector= input_vector./vector_magnitude;

%plot the vectorEC in GREEN
quiver(pointE(1,1), pointE(1,2), vectorEC(1,1), vectorEC(1,2), 'Color', 'g', 'LineWidth', 4);

%plot the vectorES in BLUE
quiver(pointE(1,1), pointE(1,2), vectorES(1,1), vectorES(1,2), 'Color', 'b', 'LineWidth', 4);

% Calculate unit orthogonal vector of vectorES
unit_orthogonal_vectorES = fcn_geometry_calcOrthogonalVector(vectorES);

% Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
quiver(pointE(1,1), pointE(1,2), unit_orthogonal_vectorES(1,1), unit_orthogonal_vectorES(1,2), 'Color', 'r', 'LineWidth', 3);

% Find the distance between the arc center and start point of the line
% segment by calculating the dot product of VectorEC and unit_orthogonal_vectorES
dist_btw_pointC_and_pointS = dot(vectorEC, unit_orthogonal_vectorES);

% Radius of the regression  Arc
radiusArc = endPointsCell{1}.fitParameters(1,3);

if dist_btw_pointC_and_pointS > radiusArc
    disp('No Intersection')
end

% One intersection point, if 
if dist_btw_pointC_and_pointS == radiusArc
    
end










