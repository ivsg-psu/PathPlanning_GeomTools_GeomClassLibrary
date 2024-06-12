% script_test_fcn_geometry_findAgreementsOfPointsToCubicPoly
% Exercises the function: fcn_geometry_findAgreementsOfPointsToCubicPoly
% Revision history:
% 2024_05_29 - Aneesh Batchu
% -- wrote the code
% 2024_06_05 - Aneesh Batcu
% -- Modified the assertions for better demonstration

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

%% Test 1: a basic linear polynomial (y = 0.5x)

rng(123)

fig_num = 1121;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0.5; % Coefficient for x
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


current_combo = [1, 2, 3, 4]; 

% Select the first four points as the source points
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 0.1; 


[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));


% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize', 20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(isequal(agreement_indices,[1; 2; 3; 4]));
assert(isequal(size(agreement_indices),[4 1]));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test 2: a basic linear polynomial (y = 0.5x)

rng(123)

fig_num = 1131;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0.5; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 7; % Number of test points to generate
sigma = 0.15; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;


current_combo = [1, 11, 18, 29]; 

% Select the first four points as the source points
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1;

% Station tolerance
station_tolerance = 1; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));


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

rng(123)

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
transverse_tolerance = 0.1;

% Base point index
current_combo = []; 

% Station tolerance
station_tolerance = []; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test 2: a basic quadratic polynomial (y = 0.2x^2) (station_tolerance)

rng(123)

fig_num = 1311;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0.2; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select any four points as the source points
current_combo = [1, 3, 5, 9];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 1; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test3: a basic quadratic polynomial (y = 0.2x^2) (station_tolerance)

rng(123)

fig_num = 1311;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0.2; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select any four points as the source points
current_combo = [3, 5, 9, 12];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.3;

% Base point index
% base_point_index = 1; 

% Station tolerance
station_tolerance = 5; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

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

rng(123)

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
current_combo = [1, 2, 3, 4];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% current_combo
current_combo = [];

% Station tolerance
station_tolerance = []; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test 2: a basic cubic polynomial (y = 0.2x^3)

rng(123)

fig_num = 3111;
figure(fig_num); clf;

a = 0.2; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
current_combo = [1, 3, 7, 9];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.5;

% Base point index
% current_combo = []; 

% Station tolerance
station_tolerance = []; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test 3: a basic cubic polynomial (y = x^3) (station_tolerance)

rng(123)

fig_num = 3111;
figure(fig_num); clf;

a = 0.2; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
current_combo = [9, 12, 13, 15];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1;

% Station tolerance
station_tolerance = 0.2; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

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

rng(123)

fig_num = 1112;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
current_combo = [2, 5, 13, 18];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1;

% Base point index
current_combo = []; 

% Station tolerance
station_tolerance = []; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));


%% Test 2: a basic cubic polynomial (y = 5)

rng(123)

fig_num = 1113;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 4; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
current_combo = [1, 3, 7, 9];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.2;

% Base point index
current_combo = []; 

% Station tolerance
station_tolerance = []; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% Test 3: a basic cubic polynomial (y = 5) (station_tolerance)

rng(123)

fig_num = 1114;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 5; % Constant term
x_range = [-2, 2]; % Range of x values
M = 4; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 1;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Did not corrupt the test points to get the same result everytime
points = corrupted_test_points;

% Select first four points as the source points
current_combo = [5, 8, 11, 15];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1; 

% Station tolerance
station_tolerance = 0.1; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

%% MISC - 1

rng(123)

fig_num = 22221;
figure(fig_num); clf;


a = 0.01; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-6, -0.5]; % Range of x values
M = 4; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points1 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


a = 0; % Coefficient for x^3
b = -0.1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 2; % Constant term
x_range = [0, 4]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points2 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


points = [corrupted_test_points1; corrupted_test_points2]; 

% Select first four points as the source points
current_combo = [2, 8, 10, 43];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1; 

% Station tolerance
station_tolerance = 0.1;

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));


%% MISC - 2

rng(123)

fig_num = 22221;
figure(fig_num); clf;


a = 0.01; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-6, -0.5]; % Range of x values
M = 4; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points1 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


a = 0; % Coefficient for x^3
b = -0.1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 2; % Constant term
x_range = [4, 8]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points2 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


points = [corrupted_test_points1; corrupted_test_points2]; 

% Select first four points as the source points
current_combo = [14, 17, 21, 26];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1; 

% Station tolerance
station_tolerance = 0.1; 

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (current_combo), (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));


%% Working - shows agreement without transverse agreement check


rng(123)

fig_num = 22221;
figure(fig_num); clf;


a = 0.01; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-6, -0.5]; % Range of x values
M = 1; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points1 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


a = 0; % Coefficient for x^3
b = -0.1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 2; % Constant term
x_range = [0, 4]; % Range of x values
M = 1; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points2 = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


points = [corrupted_test_points1; corrupted_test_points2]; 

% Select first four points as the source points
current_combo = [1, 2, 4, 11];
test_source_points = points(current_combo,:);

% Find the fitted parameters using polyfit
fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

% Tranverse tolerance
transverse_tolerance = 0.1; 

% Station tolerance
station_tolerance = [];

[agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (station_tolerance), (fig_num));

% plot the source points of the cubic polynomial curve
plot(test_source_points(:,1), test_source_points(:,2), 'b.', 'MarkerSize',20)

% Plot the fitted polynomial
x_fit = linspace(min(test_source_points(:,1)), max(test_source_points(:,1)), 100);
y_fit = polyval(fittedParameters, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

assert(length(agreement_indices(:,1))>1);
assert(isequal(size(agreement_indices,2),1));
assert(length(dist_btw_points_and_cubic_curve(:,1))>1);
assert(isequal(size(dist_btw_points_and_cubic_curve,2),1));

