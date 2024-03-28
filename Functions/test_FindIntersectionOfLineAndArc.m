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

%% In this case, one arc and one line have been used as the test data. (No Intersection)

% % Fill data points with lines and arcs
% rng(3423)
% 
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% % fig_num = 111;
% fig_num = -1;
% 
% % Arc 1 test points
% % seed_points = [1 1; 2 1.8; 3 3];
% seed_points = [2 3; 4 5; 6 3];
% M = 10; % Points per meter
% sigma = 0.02;
% 
% % Fill test data for arc 1
% [test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
% 
% fig_num = 112;
% corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 
% % Line 1 test points
% seed_points = [5 6; 9 3];
% M = 10;
% sigma = 0.02;
% 
% % fig_num = 111;
% fig_num = -1;
% % Fill test data for line 1
% test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);
% 
% 
% fig_num = 112;
% corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);

%%  In this case, one arc and one line have been used as the test data. (Two Intersection points)

% % Fill data points with lines and arcs
% rng(3423)
% 
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% % fig_num = 111;
% fig_num = -1;
% 
% % Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
% % seed_points = [4 0; 0 4; -4 0];
% M = 10; % Points per meter
% sigma = 0.02;
% 
% % Fill test data for arc 1
% [test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
% 
% fig_num = 112;
% corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 
% % Line 1 test points
% seed_points = [5 4; 8 3];
% M = 10;
% sigma = 0.02;
% 
% % fig_num = 111;
% fig_num = -1;
% % Fill test data for line 1
% test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);
% 
% 
% fig_num = 112;
% corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);

%% In this case, one arc and one line have been used as the test data. (Two Intersection points)

% % Fill data points with lines and arcs
% rng(3423)
% 
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% % fig_num = 111;
% fig_num = -1;
% 
% % Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
% % seed_points = [2 3; 4 5; 6 3];
% % seed_points = [4 0; 0 4; -4 0];
% M = 10; % Points per meter
% sigma = 0.02;
% 
% % Fill test data for arc 1
% [test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
% 
% fig_num = 112;
% corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 
% % Line 1 test points
% seed_points = [7 3; 9 4];
% M = 10;
% sigma = 0.02;
% 
% % fig_num = 111;
% fig_num = -1;
% % Fill test data for line 1
% test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);
% 
% 
% fig_num = 112;
% corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);


%% In this case, one arc and one line have been used as the test data. (One Intersection point)

% % Fill data points with lines and arcs
% rng(3423)
% 
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% % fig_num = 111;
% fig_num = -1;
% 
% % Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
% % seed_points = [2 3; 4 5; 6 3];
% % seed_points = [4 0; 0 4; -4 0];
% M = 10; % Points per meter
% sigma = 0.02;
% 
% % Fill test data for arc 1
% [test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
% 
% fig_num = 112;
% corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 
% % Line 1 test points
% seed_points = [7 4; 9 4];
% M = 10;
% sigma = 0.02;
% 
% % fig_num = 111;
% fig_num = -1;
% % Fill test data for line 1
% test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);
% 
% 
% fig_num = 112;
% corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 

% %% In this case, one arc and one line have been used as the test data. (Two Intersection points)
% 
% % Fill data points with lines and arcs
% rng(3423)
% 
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% % fig_num = 111;
% fig_num = -1;
% 
% % Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
% % seed_points = [2 3; 4 5; 6 3];
% % seed_points = [4 0; 0 4; -4 0];
% M = 10; % Points per meter
% sigma = 0.02;
% 
% % Fill test data for arc 1
% [test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);
% 
% fig_num = 112;
% corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);
% 
% % Line 1 test points
% seed_points = [4 5; 7 5];
% M = 10;
% sigma = 0.02;
% 
% % fig_num = 111;
% fig_num = -1;
% % Fill test data for line 1
% test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);
% 
% 
% fig_num = 112;
% corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
%     (probability_of_corruption), (magnitude_of_corruption), fig_num);

%% In this case, one arc and one line have been used as the test data. (One Intersection points)

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% fig_num = 111;
fig_num = -1;

% Arc 1 test points
seed_points = [2 3; 4 5; 6 3];
% seed_points = [2 3; 4 5; 6 3];
% seed_points = [4 0; 0 4; -4 0];
M = 10; % Points per meter
sigma = 0.02;

% Fill test data for arc 1
[test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

fig_num = 112;
corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Line 1 test points
seed_points = [4 1; 7 1];
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

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, fig_num);


hold on 

% Location of S (Start Point of Line Segment) 
pointS = endPointsCell{2}.firstEndPoint;

% Location of E (End Point of Line Segement)
pointE = endPointsCell{2}.lastEndPoint;

% Location of Arc center
pointC = endPointsCell{1}.fitParameters(1,1:2);

% Plot Arc Center
plot(pointC(1,1), pointC(1,2), 'r.','MarkerSize',30)

% Calculate vector EC
vectorEC = pointC - pointE;

% Calculate vector ES
vectorES = pointS - pointE;

% vector_magnitude = sum(input_vector.^2,2).^0.5;
% unitVector_input_vector= input_vector./vector_magnitude;

% Plot the vectorEC in GREEN
quiver(pointE(1,1), pointE(1,2), vectorEC(1,1), vectorEC(1,2), 'Color', 'g', 'LineWidth', 4);

% Plot the vectorES in BLUE
quiver(pointE(1,1), pointE(1,2), vectorES(1,1), vectorES(1,2), 'Color', 'b', 'LineWidth', 4);

% Calculate the unit vector of vectorES
vectorES_magnitude = sum(vectorES.^2,2).^0.5;
unit_vectorES = vectorES./vectorES_magnitude;

% Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
quiver(pointE(1,1), pointE(1,2), unit_vectorES(1,1), unit_vectorES(1,2), 'Color', 'y', 'LineWidth', 3);

% Calculate unit orthogonal vector of vectorES
unit_orthogonal_vectorES = unit_vectorES*[0 1; -1 0];

% Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
quiver(pointE(1,1), pointE(1,2), unit_orthogonal_vectorES(1,1), unit_orthogonal_vectorES(1,2), 'Color', 'r', 'LineWidth', 3);

% Find the distance between the arc center and start point of the line
% segment by calculating the dot product of VectorEC and unit_orthogonal_vectorES
dist_btw_pointC_and_pointS = dot(vectorEC, unit_orthogonal_vectorES);

% Radius of the regression Arc
radiusArc = endPointsCell{1}.fitParameters(1,3);

if dist_btw_pointC_and_pointS > radiusArc
    disp('No Intersection')
end

tole = 0.1;
% One intersection point, if the distance between arc center and start
% point of segment is equal to radius
if dist_btw_pointC_and_pointS <= radiusArc + tole && dist_btw_pointC_and_pointS >= radiusArc - tole
    intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    disp(intersectionPoints)
    plot(intersectionPoints(1,1), intersectionPoints(1,2), '.', 'Color', 'c', 'MarkerSize',30)
end


% Two intersection points, if the distance between arc center and start
% point of segment is less than others
if dist_btw_pointC_and_pointS < radiusArc - tole

    % Center point of the two intersection points
    % centerPoint_of_intersectionPoints = vectorES + dot(vectorEC, unit_vectorES)*unit_vectorES;
    centerPoint_of_intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    plot(centerPoint_of_intersectionPoints(1,1), centerPoint_of_intersectionPoints(1,2), '.', 'Color', 'g', 'MarkerSize',30);

    % Opposite side of right angle triangle: dist_btw_pointC_and_pointS
    dist_btw_centerPoint_and_pointC = dot(vectorEC, unit_orthogonal_vectorES);

    % Adjacent side of the right angle triangle
    dist_btw_centerPoint_and_intersectionPoints = (radiusArc^2 - dist_btw_centerPoint_and_pointC^2)^0.5;
    
    %
    intersectionPoint1 = centerPoint_of_intersectionPoints + dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoint2 = centerPoint_of_intersectionPoints - dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoints = [intersectionPoint1; intersectionPoint2]; 
    disp(intersectionPoints)
    plot(intersectionPoints(:,1), intersectionPoints(:,2), '.', 'Color', 'c', 'MarkerSize',30)
end


% Verify the points

distance_btw_interesectionPoints_center = sum([(intersectionPoints(:,1) - pointC(:,1)*ones(size(intersectionPoints,1),1)).^2, (intersectionPoints(:,2) - pointC(:,2)*ones(size(intersectionPoints,1),1)).^2],2).^0.5;

verifyThis = distance_btw_interesectionPoints_center - radiusArc.*ones(size(distance_btw_interesectionPoints_center,1),1);

if (verifyThis < 0.0000001)
    disp('Intersection Points are Verified')
end

% sum([(centerPoint_of_intersectionPoints(:,1) - pointC(:,1)).^2, (centerPoint_of_intersectionPoints(:,2) - pointC(:,2)).^2],2).^0.5;


