
close all

%% In this case, one arc and one line have been used as the test data. (One Intersection point(s))

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% fig_num = 111;
fig_num = -1;

% Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
seed_points = [10 2; 12 4; 10 6];
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
seed_points = [6 2; 8.5 2];
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
station_tolerance = 1; % Units are meters. 
threshold_max_points = 20;
input_points = [corrupted_test_points_arc1; corrupted_test_points_line1];

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 601;
% fig_num = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], fig_num);
% fcn_geometry_plotFitDomains(regression_domains, fig_num+2);

fig_num = 114;
fcn_geometry_plotFitDomains(regression_domains, fig_num);
tolerance = [];

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, fig_num);

tole = 1;
intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(endPointsCell, tole, fig_num);

% Verify the points

% Location of Arc center
pointC = endPointsCell{2}.fitParameters(1,1:2);

% Radius of the regression Arc
radiusArc = endPointsCell{2}.fitParameters(1,3);

distance_btw_interesectionPoints_center = sum([(intersectionPoints(:,1) - pointC(:,1)*ones(size(intersectionPoints,1),1)).^2, (intersectionPoints(:,2) - pointC(:,2)*ones(size(intersectionPoints,1),1)).^2],2).^0.5;

verifyThis = distance_btw_interesectionPoints_center - radiusArc.*ones(size(distance_btw_interesectionPoints_center,1),1);

if (verifyThis < 0.0000001)
    disp('Intersection Points are Verified')
end

% sum([(centerPoint_of_intersectionPoints(:,1) - pointC(:,1)).^2, (centerPoint_of_intersectionPoints(:,2) - pointC(:,2)).^2],2).^0.5;



%% In this case, one arc and one line have been used as the test data. (One Intersection point(s))

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% fig_num = 221;
fig_num = -1;

% Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
seed_points = [10 2; 12 4; 10 6];
% seed_points = [2 3; 4 5; 6 3];
% seed_points = [4 0; 0 4; -4 0];
M = 10; % Points per meter
sigma = 0.02;

% Fill test data for arc 1
[test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

fig_num = 222;
corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Line 1 test points
seed_points = [6 2; 9 2];
M = 10;
sigma = 0.02;

% fig_num = 221;
fig_num = -1;
% Fill test data for line 1
test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


fig_num = 222;
corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Hough Segmentation
% fig_num = 205;
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 20;
input_points = [corrupted_test_points_arc1; corrupted_test_points_line1];

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 206;
% fig_num = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], fig_num);
% fcn_geometry_plotFitDomains(regression_domains, fig_num+2);


fig_num = 224;
fcn_geometry_plotFitDomains(regression_domains, fig_num);
tolerance = [];

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, fig_num);

tole = 1;
intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(endPointsCell, tole, fig_num);

% Verify the points

% Location of Arc center
pointC = endPointsCell{2}.fitParameters(1,1:2);

% Radius of the regression Arc
radiusArc = endPointsCell{2}.fitParameters(1,3);

distance_btw_interesectionPoints_center = sum([(intersectionPoints(:,1) - pointC(:,1)*ones(size(intersectionPoints,1),1)).^2, (intersectionPoints(:,2) - pointC(:,2)*ones(size(intersectionPoints,1),1)).^2],2).^0.5;

verifyThis = distance_btw_interesectionPoints_center - radiusArc.*ones(size(distance_btw_interesectionPoints_center,1),1);

if (verifyThis < 0.0000001)
    disp('Intersection Points are Verified')
end


%% In this case, one arc and one line have been used as the test data. (One Intersection point(s))

% Fill data points with lines and arcs
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% fig_num = 321;
fig_num = -1;

% Arc 1 test points
% seed_points = [2 3; 4 5; 6 3];
seed_points = [6 2; 4 4; 6 6];
% seed_points = [2 3; 4 5; 6 3];
% seed_points = [4 0; 0 4; -4 0];
M = 10; % Points per meter
sigma = 0.02;

% Fill test data for arc 1
[test_points_arc1, true_circle_centers_arc1, true_circle_radii_arc1] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, fig_num);

fig_num = 322;
corrupted_test_points_arc1 = fcn_geometry_corruptPointsWithOutliers(test_points_arc1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Line 1 test points
seed_points = [6 2; 9 2];
M = 10;
sigma = 0.02;

% fig_num = 321;
fig_num = -1;
% Fill test data for line 1
test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


fig_num = 322;
corrupted_test_points_line1 = fcn_geometry_corruptPointsWithOutliers(test_points_line1,...
    (probability_of_corruption), (magnitude_of_corruption), fig_num);

% Hough Segmentation
% fig_num = 305;
fig_num = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 20;
input_points = [corrupted_test_points_arc1; corrupted_test_points_line1];

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 306;
% fig_num = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], fig_num);
% fcn_geometry_plotFitDomains(regression_domains, fig_num+2);


fig_num = 324;
fcn_geometry_plotFitDomains(regression_domains, fig_num);
tolerance = [];

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, fig_num);

tole = 1;
intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(endPointsCell, tole, fig_num);

% Verify the points

% Location of Arc center
pointC = endPointsCell{2}.fitParameters(1,1:2);

% Radius of the regression Arc
radiusArc = endPointsCell{2}.fitParameters(1,3);

distance_btw_interesectionPoints_center = sum([(intersectionPoints(:,1) - pointC(:,1)*ones(size(intersectionPoints,1),1)).^2, (intersectionPoints(:,2) - pointC(:,2)*ones(size(intersectionPoints,1),1)).^2],2).^0.5;

verifyThis = distance_btw_interesectionPoints_center - radiusArc.*ones(size(distance_btw_interesectionPoints_center,1),1);

if (verifyThis < 0.0000001)
    disp('Intersection Points are Verified')
end

