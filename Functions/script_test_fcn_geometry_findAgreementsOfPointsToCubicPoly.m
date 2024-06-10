% script_test_fcn_geometry_findAgreementsOfPointsToCubicPoly
% Exercises the function: fcn_geometry_findAgreementsOfPointsToCubicPoly
% Revision history:
% 2024_05_29 - Aneesh Batchu
% -- wrote the code

close all

%% Assertions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             _   _
%     /\                     | | (_)
%    /  \   ___ ___  ___ _ __| |_ _  ___  _ __  ___
%   / /\ \ / __/ __|/ _ \ '__| __| |/ _ \| '_ \/ __|
%  / ____ \\__ \__ \  __/ |  | |_| | (_) | | | \__ \
% /_/    \_\___/___/\___|_|   \__|_|\___/|_| |_|___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Assertions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test linear polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _      _                         _____      _                             _       _
% | |    (_)                       |  __ \    | |                           (_)     | |
% | |     _ _ __   ___  __ _ _ __  | |__) |__ | |_   _ _ __   ___  _ __ ___  _  __ _| |
% | |    | | '_ \ / _ \/ _` | '__| |  ___/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
% | |____| | | | |  __/ (_| | |    | |  | (_) | | |_| | | | | (_) | | | | | | | (_| | |
% |______|_|_| |_|\___|\__,_|_|    |_|   \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
%                                                 __/ |
%                                                |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Linear%20Polynomial%20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 1: a basic linear polynomial (y = 2x)

fig_num = 1121;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 2; % Coefficient for x
d = 0; % Constant term
x_range = [-1, 1]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select the first four points as the source points
test_source_points = points(1:4,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 0.1; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,[1; 2; 3]));
assert(isequal(size(agreement_indices),[3 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 2: a basic linear polynomial (y = 2x)

fig_num = 1131;
figure(fig_num); clf;

points_required_for_agreement = 10;
 
flag_find_only_best_agreement = 0; 
flag_do_plots = 1;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 2; % Coefficient for x
d = 0; % Constant term
x_range = [-1, 1]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;


% Select any four points as the source points
combo = [1, 3, 4, 5];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 0.1; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,[1; 2; 3; 4]));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test Quadratic Polynomial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ____                  _           _   _        _____      _                             _       _
%  / __ \                | |         | | (_)      |  __ \    | |                           (_)     | |
% | |  | |_   _  __ _  __| |_ __ __ _| |_ _  ___  | |__) |__ | |_   _ _ __   ___  _ __ ___  _  __ _| |
% | |  | | | | |/ _` |/ _` | '__/ _` | __| |/ __| |  ___/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
% | |__| | |_| | (_| | (_| | | | (_| | |_| | (__  | |  | (_) | | |_| | | | | (_) | | | | | | | (_| | |
%  \___\_\\__,_|\__,_|\__,_|_|  \__,_|\__|_|\___| |_|   \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
%                                                                __/ |
%                                                               |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Quadratic%20Polynomial%20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 1: a basic quadratic polynomial (y = x^2)

fig_num = 1211;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select the first four points as the source points
test_source_points = points(1:4,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = []; 

% Station tolerance
station_tolerance = []; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,[1; 2; 3; 4]));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 2: a basic quadratic polynomial (y = x^2) (station_tolerance)

fig_num = 1311;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select any four points as the source points
combo = [1, 3, 5, 9];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 0.5; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(1:8)'));
assert(isequal(size(agreement_indices),[8 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test test: a basic quadratic polynomial (y = x^2) (station_tolerance)

rng(123)

fig_num = 1311;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 3; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select any four points as the source points
combo = [3, 5, 9, 12];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.3;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 5; 

total_points_including_source_points = 20;

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(3:11)'));
assert(isequal(size(agreement_indices),[9 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 3: a basic quadratic polynomial (y = x^2) (station_tolerance)

fig_num = 1411;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select any four points as the source points
combo = [3, 5, 6, 7];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 3; 

% Station tolerance
station_tolerance = 0.5; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(3:6)'));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test Cubic Polynomial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____      _     _        _____      _                             _       _
%  / ____|    | |   (_)      |  __ \    | |                           (_)     | |
% | |    _   _| |__  _  ___  | |__) |__ | |_   _ _ __   ___  _ __ ___  _  __ _| |
% | |   | | | | '_ \| |/ __| |  ___/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
% | |___| |_| | |_) | | (__  | |  | (_) | | |_| | | | | (_) | | | | | | | (_| | |
%  \_____\__,_|_.__/|_|\___| |_|   \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
%                                           __/ |
%                                          |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Cubic%20Polynomial%20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 1: a basic cubic polynomial (y = x^3)

fig_num = 2111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [1, 2, 3, 4];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = []; 

% Station tolerance
station_tolerance = []; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,[1; 2; 3; 4]));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 2: a basic cubic polynomial (y = x^3)

fig_num = 3111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [1, 3, 7, 9];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = []; 

% Station tolerance
station_tolerance = []; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(1:9)'));
assert(isequal(size(agreement_indices),[9 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 2: a basic cubic polynomial (y = x^3) (station_tolerance)

fig_num = 3111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2.1, 2.1]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [1, 3, 7, 9];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 1.4; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(1:8)'));
assert(isequal(size(agreement_indices),[8 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test constant polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                _              _     _____      _                             _       _
%  / ____|              | |            | |   |  __ \    | |                           (_)     | |
% | |     ___  _ __  ___| |_ __ _ _ __ | |_  | |__) |__ | |_   _ _ __   ___  _ __ ___  _  __ _| |
% | |    / _ \| '_ \/ __| __/ _` | '_ \| __| |  ___/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
% | |___| (_) | | | \__ \ || (_| | | | | |_  | |  | (_) | | |_| | | | | (_) | | | | | | | (_| | |
%  \_____\___/|_| |_|___/\__\__,_|_| |_|\__| |_|   \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|
%                                                           __/ |
%                                                          |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Linear%20Polynomial%20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 1: a basic cubic polynomial (y = 5)

fig_num = 1112;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [1, 2, 3, 4];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.2;

% Base point index
% base_point_index = []; 

% Station tolerance
station_tolerance = []; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,[1; 2; 3; 4]));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 2: a basic cubic polynomial (y = 5)

fig_num = 1113;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [1, 3, 7, 9];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.2;

% Base point index
% base_point_index = []; 

% Station tolerance
station_tolerance = []; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(1:9)'));
assert(isequal(size(agreement_indices),[9 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));

%% Test 3: a basic cubic polynomial (y = 5) (station_tolerance)

fig_num = 1114;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 0;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
combo = [2, 3, 5, 7];
test_source_points = points(combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.2;

% Base point index
% base_point_index = 3; 

% Station tolerance
station_tolerance = 0.1; 

total_points_including_source_points = 20; 

% [agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (base_point_index), (station_tolerance), (fig_num));
[agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num));

assert(isequal(agreement_indices,(2:6)'));
assert(isequal(size(agreement_indices),[5 1]));
assert(length(polygon_vertices(:,1))>1);
assert(isequal(size(polygon_vertices,2),2));
