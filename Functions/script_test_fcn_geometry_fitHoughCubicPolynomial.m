% script_test_fcn_geometry_fitHoughCubicPolynomial
% Exercises the function: fcn_geometry_fitHoughCubicPolynomial
% Revision history:
% 2024_05_16 - Aneesh Batchu
% -- wrote the code
% 2024_05_20 - Aneesh Batchu
% -- Added assertions


% To-do
%
% 2024_05_20 - Aneesh Batchu
% -- Commented the Test: 4 a basic constant (y = 5) because of warnings.
% The warnings are displayed due to
% "fcn_geometry_findAgreementOfPointsToCubicPoly". Change the fig_num = -1
% to some number (123) to show the agreement plots. 

close all;

%% Test 1: a basic cubic polynomial (y = x^3)

rng(123)

% fig_num = 111;
% figure(fig_num);
% clf;

fig_num = -1;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.15; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

% fig_num = 1111;
% figure(fig_num); 
% clf;
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11111;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = 5;

flag_find_only_best_agreement = 1; 

station_tolerance = 0.1;

domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

% Check the regression fit
% regression_domains = fcn_geometry_HoughRegression(domains, transverse_tolerance, fig_num);

% Check the output type and size
for ith_domain = 1:length(domains)-1
    domain = domains{ith_domain};
    assert(isstruct(domain));
    assert(isfield(domain,'best_fit_type'));
    assert(isfield(domain,'points_in_domain'));
    assert(isfield(domain,'best_fit_parameters'));
    assert(isfield(domain,'best_fit_domain_box'));
    assert(isfield(domain,'best_fit_1_sigma_box'));
    assert(isfield(domain,'best_fit_2_sigma_box'));
    assert(isfield(domain,'best_fit_3_sigma_box'));
    assert(ischar(domain.best_fit_type));
    assert(strcmp('Hough cubic polynomial',domain.best_fit_type));
    assert(length(domain.points_in_domain(:,1))>1);
    assert(length(domain.points_in_domain(1,:))==2);
    assert(isequal(size(domain.best_fit_parameters),[1 4]));
end

% Check the last domain (unfitted points)
domain = domains{end};
assert(isstruct(domain));
assert(isfield(domain,'best_fit_type'));
assert(isfield(domain,'points_in_domain'));
assert(isfield(domain,'best_fit_parameters'));
assert(isfield(domain,'best_fit_domain_box'));
assert(ischar(domain.best_fit_type));
assert(strcmp('unfitted',domain.best_fit_type));
assert(length(domain.points_in_domain(:,1))>1);
assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));

%% Test 2: a basic quadratic polynomial (y = x^2)

rng(123)

% fig_num = 112;
% figure(fig_num);
% clf;

fig_num = -1;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 0]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.8; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

% fig_num = 1112;
% figure(fig_num); 
% clf;
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11112;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = [];

flag_find_only_best_agreement = []; 

station_tolerance = 0.5;

domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

% Check the output type and size
for ith_domain = 1:length(domains)-1
    domain = domains{ith_domain};
    assert(isstruct(domain));
    assert(isfield(domain,'best_fit_type'));
    assert(isfield(domain,'points_in_domain'));
    assert(isfield(domain,'best_fit_parameters'));
    assert(isfield(domain,'best_fit_domain_box'));
    assert(isfield(domain,'best_fit_1_sigma_box'));
    assert(isfield(domain,'best_fit_2_sigma_box'));
    assert(isfield(domain,'best_fit_3_sigma_box'));
    assert(ischar(domain.best_fit_type));
    assert(strcmp('Hough cubic polynomial',domain.best_fit_type));
    assert(length(domain.points_in_domain(:,1))>1);
    assert(length(domain.points_in_domain(1,:))==2);
    assert(isequal(size(domain.best_fit_parameters),[1 4]));
end

% Check the last domain (unfitted points)
domain = domains{end};
assert(isstruct(domain));
assert(isfield(domain,'best_fit_type'));
assert(isfield(domain,'points_in_domain'));
assert(isfield(domain,'best_fit_parameters'));
assert(isfield(domain,'best_fit_domain_box'));
assert(ischar(domain.best_fit_type));
assert(strcmp('unfitted',domain.best_fit_type));
assert(length(domain.points_in_domain(:,1))>1);
assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));

%% Test 3: a basic linear polynomial (y = x)

rng(123)

% fig_num = 113;
% figure(fig_num);
% clf;

fig_num = -1;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 1; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

% fig_num = 1113;
% figure(fig_num); 
% clf;
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11113;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.1;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 

station_tolerance = 0.05;

domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

% Check the output type and size
for ith_domain = 1:length(domains)-1
    domain = domains{ith_domain};
    assert(isstruct(domain));
    assert(isfield(domain,'best_fit_type'));
    assert(isfield(domain,'points_in_domain'));
    assert(isfield(domain,'best_fit_parameters'));
    assert(isfield(domain,'best_fit_domain_box'));
    assert(isfield(domain,'best_fit_1_sigma_box'));
    assert(isfield(domain,'best_fit_2_sigma_box'));
    assert(isfield(domain,'best_fit_3_sigma_box'));
    assert(ischar(domain.best_fit_type));
    assert(strcmp('Hough cubic polynomial',domain.best_fit_type));
    assert(length(domain.points_in_domain(:,1))>1);
    assert(length(domain.points_in_domain(1,:))==2);
    assert(isequal(size(domain.best_fit_parameters),[1 4]));
end

% Check the last domain (unfitted points)
domain = domains{end};
assert(isstruct(domain));
assert(isfield(domain,'best_fit_type'));
assert(isfield(domain,'points_in_domain'));
assert(isfield(domain,'best_fit_parameters'));
assert(isfield(domain,'best_fit_domain_box'));
assert(ischar(domain.best_fit_type));
assert(strcmp('unfitted',domain.best_fit_type));
assert(length(domain.points_in_domain(:,1))>1);
assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));

%% Test 4: a basic constant (y = 5): WARNINGS

% fig_num = 114;
% figure(fig_num);
% clf;
% 
% a = 0; % Coefficient for x^3
% b = 0; % Coefficient for x^2
% c = 0; % Coefficient for x
% d = 1; % Constant term
% x_range = [-2, 2]; % Range of x values
% M = 10; % Number of test points to generate
% sigma = 0; % Standard deviation for randomness
% 
% [test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);
% 
% fig_num = 1114;
% figure(fig_num); 
% clf;
% % Corrupt the results
% probability_of_corruption = 0.1;
% magnitude_of_corruption = 3;
% 
% corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
%     (probability_of_corruption), (magnitude_of_corruption), (fig_num));
% 
% fig_num = 11114;
% figure(fig_num); clf;
% 
% inputPoints = corrupted_test_points;
% transverse_tolerance = 0.3;
% 
% points_required_for_agreement = [];
% 
% flag_find_only_best_agreement = 1; 
% 
% station_tolerance = 0.1;
% 
% domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));
% 
% % Check the output type and size
% for ith_domain = 1:length(domains)-1
%     domain = domains{ith_domain};
%     assert(isstruct(domain));
%     assert(isfield(domain,'best_fit_type'));
%     assert(isfield(domain,'points_in_domain'));
%     assert(isfield(domain,'best_fit_parameters'));
%     assert(isfield(domain,'best_fit_domain_box'));
%     assert(isfield(domain,'best_fit_1_sigma_box'));
%     assert(isfield(domain,'best_fit_2_sigma_box'));
%     assert(isfield(domain,'best_fit_3_sigma_box'));
%     assert(ischar(domain.best_fit_type));
%     assert(strcmp('Hough cubic polynomial',domain.best_fit_type));
%     assert(length(domain.points_in_domain(:,1))>1);
%     assert(length(domain.points_in_domain(1,:))==2);
%     assert(isequal(size(domain.best_fit_parameters),[1 4]));
% end
% 
% % Check the last domain (unfitted points)
% domain = domains{end};
% assert(isstruct(domain));
% assert(isfield(domain,'best_fit_type'));
% assert(isfield(domain,'points_in_domain'));
% assert(isfield(domain,'best_fit_parameters'));
% assert(isfield(domain,'best_fit_domain_box'));
% assert(ischar(domain.best_fit_type));
% assert(strcmp('unfitted',domain.best_fit_type));
% assert(length(domain.points_in_domain(:,1))>1);
% assert(length(domain.points_in_domain(1,:))==2);
% assert(isnan(domain.best_fit_parameters));

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
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 11115;
figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 1;

points_required_for_agreement = [];

flag_find_only_best_agreement = []; 

station_tolerance = [];

domains = fcn_geometry_fitHoughCubicPolynomial(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (fig_num));

% Check the output type and size
for ith_domain = 1:length(domains)-1
    domain = domains{ith_domain};
    assert(isstruct(domain));
    assert(isfield(domain,'best_fit_type'));
    assert(isfield(domain,'points_in_domain'));
    assert(isfield(domain,'best_fit_parameters'));
    assert(isfield(domain,'best_fit_domain_box'));
    assert(isfield(domain,'best_fit_1_sigma_box'));
    assert(isfield(domain,'best_fit_2_sigma_box'));
    assert(isfield(domain,'best_fit_3_sigma_box'));
    assert(ischar(domain.best_fit_type));
    assert(strcmp('Hough cubic polynomial',domain.best_fit_type));
    assert(length(domain.points_in_domain(:,1))>1);
    assert(length(domain.points_in_domain(1,:))==2);
    assert(isequal(size(domain.best_fit_parameters),[1 4]));
end

% Check the last domain (unfitted points)
domain = domains{end};
assert(isstruct(domain));
assert(isfield(domain,'best_fit_type'));
assert(isfield(domain,'points_in_domain'));
assert(isfield(domain,'best_fit_parameters'));
assert(isfield(domain,'best_fit_domain_box'));
assert(ischar(domain.best_fit_type));
assert(strcmp('unfitted',domain.best_fit_type));
assert(length(domain.points_in_domain(:,1))>1);
assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));