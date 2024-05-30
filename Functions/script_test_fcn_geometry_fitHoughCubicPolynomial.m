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

rng(123)

fig_num = 1121;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 2; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;


% Tranverse tolerance
transverse_tolerance = 0.5;

% Station tolerance
station_tolerance = [];

points_required_for_agreement = 10;

flag_find_only_best_agreement = 0; 

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
% assert(length(domain.points_in_domain(:,1))>1);
% assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));

%% Test 2: a basic linear polynomial (y = 2x) (Station Tolerance)

rng(123)

fig_num = 1131;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 2; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.1; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;


% Tranverse tolerance
transverse_tolerance = 0.2;

% Station tolerance
station_tolerance = 0.1;

points_required_for_agreement = 4;

flag_find_only_best_agreement = 0; 

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
% assert(length(domain.points_in_domain(:,1))>1);
% assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));

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


%% Test 1: a basic quadratic polynomial (y = x^2) (Station Tolerance)

rng(123)

fig_num = 1211;
figure(fig_num); clf;


a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 3; % Number of test points to generate
sigma = 0.8; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = [];

flag_find_only_best_agreement = []; 

station_tolerance = 4;

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

%% Test 2: a basic quadratic polynomial (y = x^2) (Station Tolerance)

rng(123)

fig_num = 1311;
figure(fig_num); clf;


a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 0]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.8; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;
transverse_tolerance = 0.4;

points_required_for_agreement = 10;

flag_find_only_best_agreement = 1; 

station_tolerance = 4;

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

%% Test 3: a basic quadratic polynomial (y = x^2) 

rng(123)

fig_num = 1411;
figure(fig_num); clf;


a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0.8; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;
transverse_tolerance = 0.4;

points_required_for_agreement = 4;

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

%% Test 4: a basic quadratic polynomial (y = x^2) 

rng(123)

fig_num = 1522;
figure(fig_num); clf;


a = 0; % Coefficient for x^3
b = 0.5; % Coefficient for x^2
c = 2; % Coefficient for x
d = 1; % Constant term
x_range = [0, 5]; % Range of x values
M = 3; % Number of test points to generate
sigma = 0.8; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

inputPoints = corrupted_test_points;
transverse_tolerance = 0.4;

points_required_for_agreement = 10;

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

%% Test 1: a basic cubic polynomial (y = x^3) (Station Tolerance)

rng(123)

fig_num = 2111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-3, 3]; % Range of x values
M = 2; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% fig_num = 1111;
% figure(fig_num); 
% clf;
% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% fig_num = 11111;
% figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = 4;

flag_find_only_best_agreement = 0; 

station_tolerance = 5;

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

%% Test 2: a basic cubic polynomial (y = x^3) 

rng(123)

fig_num = 3111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-1, 1]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% fig_num = 1111;
% figure(fig_num); 
% clf;
% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% fig_num = 11111;
% figure(fig_num); clf;

inputPoints = corrupted_test_points;
transverse_tolerance = 0.2;

points_required_for_agreement = 10;

flag_find_only_best_agreement = 0; 

station_tolerance = [];

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

%% Test 3: a basic cubic polynomial (y = x^3) 

rng(123)

fig_num = 4111;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-2, 2]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.2;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


inputPoints = corrupted_test_points;
transverse_tolerance = 0.2;

points_required_for_agreement = 10;

flag_find_only_best_agreement = 0; 

station_tolerance = [];

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

%% Test 1: a basic constant (y = 5) (Station tolerance)

fig_num = 1112;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 1; % Constant term
x_range = [-2, 2]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.01; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = [];

flag_find_only_best_agreement = 1; 

station_tolerance = 0.1;

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

%% Test 1: a basic constant (y = 5) 

fig_num = 1113;
figure(fig_num); clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 1; % Constant term
x_range = [-1, 1]; % Range of x values
M = 10; % Number of test points to generate
sigma = 0.01; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, -1);

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


inputPoints = corrupted_test_points;
transverse_tolerance = 0.3;

points_required_for_agreement = 5;

flag_find_only_best_agreement = 1; 

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

%% Test miscellaneous polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  __  __ _              _ _
% |  \/  (_)            | | |
% | \  / |_ ___  ___ ___| | | __ _ _ __   ___  ___  _   _ ___
% | |\/| | / __|/ __/ _ \ | |/ _` | '_ \ / _ \/ _ \| | | / __|
% | |  | | \__ \ (_|  __/ | | (_| | | | |  __/ (_) | |_| \__ \
% |_|  |_|_|___/\___\___|_|_|\__,_|_| |_|\___|\___/ \__,_|___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 5: a cubic polynomial (y = 2x^3 - x^2 + 5x + 2)

fig_num = 2224;
figure(fig_num); clf;

a = 1; % Coefficient for x^3
b = -1; % Coefficient for x^2
c = 0.2; % Coefficient for x
d = 2; % Constant term
x_range = [-4, 4]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, ~] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, (-1));

% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 2;

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));


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
% assert(length(domain.points_in_domain(:,1))>1);
% assert(length(domain.points_in_domain(1,:))==2);
assert(isnan(domain.best_fit_parameters));
