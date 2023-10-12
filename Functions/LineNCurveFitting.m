%% Line Fitting and Curve Fitting


%% Fitting a Line - by generating a random data

rng("default");

% Generate random x-values 
x = 10 * rand(100, 1);

% Create corresponding y-values with some noise
noise = 0.4 * randn(100, 1); 
slope = 2; 
intercept = 5; % Y-intercept
y = slope * x + intercept + noise;

% Fit a line (y = mx + b) to the data
p = polyfit(x, y, 1);

figure(1)
% Plot the data points
scatter(x, y, 'b', 'filled'); 

hold on; 
grid on;

% Plot the fitted line
xFit = linspace(min(x), max(x), 100); % Generate x-values for the line
yFit = polyval(p, xFit); % Compute corresponding y-values for the line
plot(xFit, yFit, 'r', 'LineWidth', 2); % Plot the fitted line


xlabel('X-axis');
ylabel('Y-axis');
title('Random Data Fitted to a Line');
legend('Data Points', 'Fitted Line', 'Location', 'best');

%% Fitting an Arc - by generating random data

rng("default");

startAngle = 0;  % Starting angle
endAngle = pi/3; % Ending angle 

% Generate random angles between startAngle and endAngle
angles = startAngle + (endAngle - startAngle)  * rand(100, 1);

% Generate random radii between 5 and 7
radii = 20 + 5 * rand(100, 1);

% Calculate x and y coordinates of the data points
x = radii .* cos(angles);  % 5 + radii .* cos(angles);
y = radii .* sin(angles);  % 3 + radii .* sin(angles);

% Add noise to the data
noise = 0.1 * randn(100, 1); 
x = x + noise;
y = y + noise;

% Define the function to fit the data to an arc
arcFit = @(parameters, angles) [parameters(1) + parameters(3) * cos(angles), parameters(2) + parameters(3) * sin(angles)];

% Initial guess for parameters: [x_center, y_center, radius]
initialGuess = [1, 1, 1];

% Fit the data to a circle
fittedParameters = lsqcurvefit(arcFit, initialGuess, angles, [x, y]);

% Extract the fitted parameters
x_center = fittedParameters(1);
y_center = fittedParameters(2);
radius = fittedParameters(3);

figure(12)
% Plot the data points
scatter(x, y, 'b', 'filled');

hold on; 
grid on

% Plot the fitted arc
fittedAngles = linspace(startAngle, endAngle, 100);
fittedCircle = arcFit(fittedParameters, fittedAngles');
plot(fittedCircle(:,1), fittedCircle(:,2), 'r', 'LineWidth', 2);

% Mark the center of the fitted arc
plot(x_center, y_center, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');


xlabel('X-axis');
ylabel('Y-axis');
title('Random Data Fitted to an Arc');
legend('Data Points', 'Fitted Arc', 'Location','best');


%% Fitting a spiral - by generating random data

rng("default");

% Generate random data points
x = linspace(10, 20, 100)';
a_true = 0.4; % True parameter for the spiral

% Generate corresponding y-values based on a spiral
y_true = a_true * x.^2;

% Add some noise to the data
noise = 0.6 * randn(size(x));
y_data = y_true + noise;

% Define the custom spiral function
spiralFunction = @(parameters, x) parameters(1) * x.^2;

% Initial guess for the parameter
initialGuess = [0.05]; 

% Fit the data to the spiral function
fitParameters = lsqcurvefit(spiralFunction, initialGuess, x, y_data);

% Plot the fitted spiral
fittedSpiral = spiralFunction(fitParameters, x);

figure(123)

% Plot the data and the fitted curve
scatter(x, y_data, 'b', 'filled');
hold on;
grid on;
plot(x, fittedSpiral, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Spiral');
xlabel('X');
ylabel('Y');
legend('Data Points', 'Fitted Spiral', 'Location','best');
title('Random Data Fitted to a Spiral');

% Display the fitted parameter
fprintf('Fitted parameter for the spiral: a = %.4f\n', fitParameters(1));

