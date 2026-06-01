% script_test_fcn_geometry_sortRegressionDomains.m
% tests fcn_geometry_sortRegressionDomains.m

% This script uses "fcn_geometry_HoughRegression.m" to fit the segments

% Revision History
% 2024_02_28 - Aneesh Batchu
% - Started the script

%% Clear workspace

close all

%% Two line segments: The ending point of a line segment is joined to the starting point of the other line segment 

rng(343)

figNum = -1;

    % Line test points
seed_points = [1 2; 2 4; 4 5];
M = 10;
sigma = 0.02;

test_points_twoLineSegments = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Corrupt the points
% figNum = 113;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(test_points_twoLineSegments,...
    (probability_of_corruption), (magnitude_of_corruption),figNum);


% Hough Segmentation
% figNum = 114;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, figNum);

clf;
figNum = 113;
fcn_geometry_plotFitDomains(domains, figNum);

% figure(figNum)

tolerance = [];

[endPointsCell, sortedHoughSegmentEndPoints, ~] = fcn_geometry_sortRegressionDomains(domains, tolerance, figNum);

endPointsCell_struct = endPointsCell{end};
assert(isstruct(endPointsCell_struct));
assert(isfield(endPointsCell_struct,'firstEndPoint'));
assert(isfield(endPointsCell_struct,'lastEndPoint'));
assert(isfield(endPointsCell_struct,'fitType'));
assert(isfield(endPointsCell_struct,'fitParameters'));
assert(ischar(endPointsCell_struct.fitType));
assert(length(sortedHoughSegmentEndPoints(:,1))>1);
assert(length(sortedHoughSegmentEndPoints(1,:))==4);


if 1==0
    disp(endPointsCell)
    disp(sortedHoughSegmentEndPoints)
   
end

% The plots would serve the purpose of assertions

%% Two line segments: The ending point of a line segment is joined to the starting point of the other line segment 

rng(343)

figNum = 112;

% Line test points
seed_points = [1 2; 2 4; 4 5];
M = 10;
sigma = 0.02;

test_points_twoLineSegments = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Corrupt the points
figNum = -1;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(test_points_twoLineSegments,...
    (probability_of_corruption), (magnitude_of_corruption),figNum);


% Hough Segmentation
% figNum = 114;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, figNum);

clf;
figNum = 113;
fcn_geometry_plotFitDomains(domains, figNum);
% figure(figNum)

tolerance = 0.2;

[endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(domains, tolerance, figNum);

endPointsCell_struct = endPointsCell{end};
assert(isstruct(endPointsCell_struct));
assert(isfield(endPointsCell_struct,'firstEndPoint'));
assert(isfield(endPointsCell_struct,'lastEndPoint'));
assert(isfield(endPointsCell_struct,'fitType'));
assert(isfield(endPointsCell_struct,'fitParameters'));
assert(ischar(endPointsCell_struct.fitType));
assert(length(sortedHoughSegmentEndPoints(:,1))>1);
assert(length(sortedHoughSegmentEndPoints(1,:))==4);
assert(length(closeEndPointsMatrix(:,1))>1);
assert(length(closeEndPointsMatrix(1,:))==2);

if 1==0
    disp(endPointsCell)
    disp(sortedHoughSegmentEndPoints)
    disp(closeEndPointsMatrix)
end

% The plots would serve the purpose of assertions
%% Two line segments: The ending point of a line segment is not joined to the starting point of the other line segment 

rng(343)

figNum = -1;

% Line 1 test points
seed_points = [1 2; 2 4];
M = 10;
sigma = 0.02;

test_points_LineSegment1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);


% Line 2 test points
seed_points = [2.5 4.5; 4 5];
M = 10;
sigma = 0.02;

test_points_LineSegment2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

testpoints = [test_points_LineSegment1; test_points_LineSegment2];

% Corrupt the points
% figNum = 116;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),figNum);


% Hough Segmentation
% figNum = 117;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 10;
input_points = corrupted_testpoints;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, figNum);

clf;
figNum = 116;
fcn_geometry_plotFitDomains(domains, figNum);
tolerance = 1;

[endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(domains, tolerance, figNum);

endPointsCell_struct = endPointsCell{end};
assert(isstruct(endPointsCell_struct));
assert(isfield(endPointsCell_struct,'firstEndPoint'));
assert(isfield(endPointsCell_struct,'lastEndPoint'));
assert(isfield(endPointsCell_struct,'fitType'));
assert(isfield(endPointsCell_struct,'fitParameters'));
assert(ischar(endPointsCell_struct.fitType));
assert(length(sortedHoughSegmentEndPoints(:,1))>1);
assert(length(sortedHoughSegmentEndPoints(1,:))==4);
assert(length(closeEndPointsMatrix(:,1))>1);
assert(length(closeEndPointsMatrix(1,:))==2);

if 1==0
    disp(endPointsCell)
    disp(sortedHoughSegmentEndPoints)
    disp(closeEndPointsMatrix)
end
% The plots would serve the purpose of assertions
%% Advanced test 1: Multiple segments were used to create a path 

rng(343)

% figNum = 901;
figNum = -1;

% Line test points - Line Segment 1
seed_points = [-15 0; -10 5];
M = 10;
sigma = 0.02;

test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 2
seed_points = [-10 6; -1 6];
M = 10;
sigma = 0.02;

test_points_line2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 3 & 4 & 5
seed_points = [-2 4; 1 10; 7 4; 13 13];
M = 10;
sigma = 0.02;

test_points_line345 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 6
seed_points = [10, 12; 18, 3];
M = 10;
sigma = 0.02;

test_points_line6 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);


testpoints = [test_points_line1; test_points_line2; test_points_line345; test_points_line6];

% Corrupt the points
% figNum = 801;
% fig_nuM = figNum;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),figNum);

% Hough Segmentation
% figNum = 701;
figNum = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 20;
input_points = corrupted_testpoints;

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, figNum);


% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], figNum);

clf;
figNum = 601;
fcn_geometry_plotFitDomains(regression_domains, figNum);

% figNum = fig_nuM;
tolerance = 50; 
% [endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(regression_domains, [], figNum);
[endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, figNum);

endPointsCell_struct = endPointsCell{end};
assert(isstruct(endPointsCell_struct));
assert(isfield(endPointsCell_struct,'firstEndPoint'));
assert(isfield(endPointsCell_struct,'lastEndPoint'));
assert(isfield(endPointsCell_struct,'fitType'));
assert(isfield(endPointsCell_struct,'fitParameters'));
assert(ischar(endPointsCell_struct.fitType));
assert(length(sortedHoughSegmentEndPoints(:,1))>1);
assert(length(sortedHoughSegmentEndPoints(1,:))==4);
assert(length(closeEndPointsMatrix(:,1))>1);
assert(length(closeEndPointsMatrix(1,:))==2);


if 1==0
    disp(endPointsCell)
    disp(sortedHoughSegmentEndPoints)
    disp(closeEndPointsMatrix)
end

% The plots would serve the purpose of assertions
%% Advanced test : Multiple segments were used to create a path 

rng(343)

% figNum = 902;
figNum = -1;

% Line test points - Line Segment 1
seed_points = [-15 0; -10 5];
M = 10;
sigma = 0.02;

test_points_line1 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 2
seed_points = [-8 6; -1 6];
M = 10;
sigma = 0.02;

test_points_line2 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 3 & 4
seed_points = [-4 8; 2 0; 8 0];
M = 10;
sigma = 0.02;

test_points_line34 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 5
seed_points = [6, -3; 10, 4];
M = 10;
sigma = 0.02;

test_points_line5 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);

% Line test points - Line Segment 6
seed_points = [12, 4; 16, 5];
M = 10;
sigma = 0.02;

test_points_line6 = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, figNum);


testpoints = [test_points_line1; test_points_line2; test_points_line34; test_points_line5; test_points_line6];

% Corrupt the points
% figNum = 802;
% fig_nuM = figNum;
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_testpoints = fcn_geometry_corruptPointsWithOutliers(testpoints,...
    (probability_of_corruption), (magnitude_of_corruption),figNum);

% Hough Segmentation
% figNum = 702;
figNum = -1;
transverse_tolerance = 0.05; % Units are meters
station_tolerance = 0.5; % Units are meters. 
threshold_max_points = 20;
input_points = corrupted_testpoints;

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, figNum);

% figNum = 602;
figNum = -1;
% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, [], figNum);

clf;
figNum = 602;
fcn_geometry_plotFitDomains(regression_domains, figNum);


% figNum = fig_nuM;
tolerance = 50; 
% [endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(regression_domains, [], figNum);
[endPointsCell, sortedHoughSegmentEndPoints, closeEndPointsMatrix] = fcn_geometry_sortRegressionDomains(regression_domains, tolerance, figNum);

endPointsCell_struct = endPointsCell{end};
assert(isstruct(endPointsCell_struct));
assert(isfield(endPointsCell_struct,'firstEndPoint'));
assert(isfield(endPointsCell_struct,'lastEndPoint'));
assert(isfield(endPointsCell_struct,'fitType'));
assert(isfield(endPointsCell_struct,'fitParameters'));
assert(ischar(endPointsCell_struct.fitType));
assert(length(sortedHoughSegmentEndPoints(:,1))>1);
assert(length(sortedHoughSegmentEndPoints(1,:))==4);
assert(length(closeEndPointsMatrix(:,1))>1);
assert(length(closeEndPointsMatrix(1,:))==2);

if 1==0

    disp(endPointsCell)
    disp(sortedHoughSegmentEndPoints)
    disp(closeEndPointsMatrix);
end
% The plots would serve the purpose of assertions
