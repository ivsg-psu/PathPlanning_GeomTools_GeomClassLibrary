% script_test_fcn_geometry_fillCubicPolyTestPoints
% Exercises the function: fcn_geometry_fillCubicPolyTestPoints
% Revision history:
% 2024_05_08 - Aneesh Batchu
% -- wrote the code

close all;

%% Test 1: a basic cubic polnomial (y = x^3)

fig_num = 1;
figure(fig_num);
clf;

a = 1; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 10; % Number of test points to generate
sigma = 2; % Standard deviation for randomness

[test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

assert(length(test_points(:,1))==length(true_points(:,1)));

% % Plot the true cubic polynomial curve and the test points
% figure;
% plot(true_points(:,1), true_points(:,2), 'b-', 'LineWidth', 2);
% hold on;
% scatter(test_points(:,1), test_points(:,2), 'r', 'filled');
% xlabel('x');
% ylabel('y');
% title('Test Points along a Cubic Polynomial Curve');
% legend('True Curve', 'Test Points');

%% Test 2: a basic quadratic polnomial (y = x^2)

fig_num = 2;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 1; % Coefficient for x^2
c = 0; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

assert(length(test_points(:,1))==length(true_points(:,1)));

%% Test 3: a basic linear polnomial (y = x)

fig_num = 3;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 1; % Coefficient for x
d = 0; % Constant term
x_range = [-5, 5]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.2; % Standard deviation for randomness

[test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

assert(length(test_points(:,1))==length(true_points(:,1)));

%% Test 4: a basic constant (y = 5)

fig_num = 4;
figure(fig_num);
clf;

a = 0; % Coefficient for x^3
b = 0; % Coefficient for x^2
c = 0; % Coefficient for x
d = 1; % Constant term
x_range = [-5, 5]; % Range of x values
M = 5; % Number of test points to generate
sigma = 0.000002; % Standard deviation for randomness

[test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

assert(length(test_points(:,1))==length(true_points(:,1)));

%% Test 5: a cubic polynomial (y = 2x^3 - x^2 + 5x + 2)

fig_num = 5;
figure(fig_num);
clf;

a = 2; % Coefficient for x^3
b = -1; % Coefficient for x^2
c = 5; % Coefficient for x
d = 2; % Constant term
x_range = [-5, 5]; % Range of x values
M = 10; % Number of test points to generate
sigma = 5; % Standard deviation for randomness

[test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints(a, b, c, d, x_range, M, sigma, fig_num);

assert(length(test_points(:,1))==length(true_points(:,1)));

