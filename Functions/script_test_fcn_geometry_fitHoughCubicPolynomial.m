% script_test_fcn_geometry_fitHoughCubicPolynomial
% Exercises the function: fcn_geometry_fitHoughCubicPolynomial
% Revision history:
% 2024_05_16 - Aneesh Batchu
% -- wrote the code

close all;

%% Test 1: a basic cubic polnomial (y = x^3)

fig_num = 111;
figure(fig_num);
clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

fig_num = 1111;
figure(fig_num); 
clf;
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 0.1;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11111;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.1;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 


domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));


domain = domains{end};
assert(isstruct(domain));

%% Test 2: a basic quadratic polnomial (y = x^2)

fig_num = 112;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 0]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

fig_num = 1112;
figure(fig_num); 
clf;
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11112;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.1;

points_required_for_agreement = [];

flag_find_only_best_agreement = []; 


domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

domain = domains{end};
assert(isstruct(domain));

%% Test 3: a basic linear polnomial (y = x)

fig_num = 113;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 1; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

fig_num = 1113;
figure(fig_num); 
clf;
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11113;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.2;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 


domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

domain = domains{end};
assert(isstruct(domain));

%% Test 4: a basic constant (y = 5)

fig_num = 114;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 1; % Constant term
x_range = [-1, 1]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

fig_num = 1114;
figure(fig_num); 
clf;
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 0.2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11114;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.2;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 


domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

domain = domains{end};
assert(isstruct(domain));

%% Test 5: a cubic polynomial (y = 2x^3 - x^2 + 5x + 2)

fig_num = 115;
figure(fig_num);
clf;

a = 2; % Coefficient for x^3
b = -1; % Coefficient for x^2
c = 5; % Coefficient for x
d = 2; % Constant term
x_range = [-5, 5]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

fig_num = 1115;
figure(fig_num); 
clf;
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 0.2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11115;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.2;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 


domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

domain = domains{end};
assert(isstruct(domain));