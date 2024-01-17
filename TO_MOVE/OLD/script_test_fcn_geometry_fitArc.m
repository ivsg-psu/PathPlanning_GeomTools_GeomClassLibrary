% script_test_fcn_geometry_fitArc.m
% tests fcn_geometry_fitArc.m

% Revision history
% 2023_10_18 - Aneesh Batchu
% -- wrote the code originally

%% Set up the workspace

clc
close all

%% Examples for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ______                           _
% |  ____|                         | |
% | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                            | |
%                            |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Arc - Case 1

% Generate synthetic noisy data resembling a circular arc
center_xCoor_true = 2;
center_yCoor_true = 3;
radius_true = 5;
startAngle_true = pi/2;  % 90 degrees in radians
endAngle_true = pi;  % 180 degrees in radians

theta_true = linspace(startAngle_true, endAngle_true, 100);
x_true = center_xCoor_true + radius_true * cos(theta_true);
x_true = x_true';
y_true = center_yCoor_true + radius_true * sin(theta_true);
y_true = y_true';

% Add some noise to the data
noise = 0.2 * randn(size(x_true,1),1);
x_noisy = x_true + noise;
y_noisy = y_true + noise;

% Initial parameter estimates (these are just initial guesses)
initial_center_xCoor_guess = 1;
initial_center_yCoor_guess = 1;
initial_radius_guess = 1;
initial_startAngle_guess = 1;
initial_endAngle_guess = 1;

% Combine noisy data into a single vector of points
noisyData = [x_noisy, y_noisy];

initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

[fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, 1);
disp(fittedArc);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Arc', 'Location', 'best');
%% Arc - Case 2

% Generate synthetic noisy data resembling a circular arc
center_xCoor_true = 0;
center_yCoor_true = 0;
radius_true = 5;
startAngle_true = pi/4;  % 45 degrees in radians
endAngle_true = pi/2;  % 90 degrees in radians

theta_true = linspace(startAngle_true, endAngle_true, 100);
x_true = center_xCoor_true + radius_true * cos(theta_true);
x_true = x_true';
y_true = center_yCoor_true + radius_true * sin(theta_true);
y_true = y_true';

% Add some noise to the data
noise = 0.2 * randn(size(x_true,1),1);
x_noisy = x_true + noise;
y_noisy = y_true + noise;

% Initial parameter estimates (these are just initial guesses)
initial_center_xCoor_guess = 5;
initial_center_yCoor_guess = 5;
initial_radius_guess = 5;
initial_startAngle_guess = 1;
initial_endAngle_guess = 1;

% Combine noisy data into a single vector of points
noisyData = [x_noisy, y_noisy];

initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

[fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, 12);
disp(fittedArc);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Arc', 'Location', 'best');

%% Arc - Case 3

% Generate synthetic noisy data resembling a circular arc
center_xCoor_true = 25;
center_yCoor_true = 15;
radius_true = 2;
startAngle_true = pi/4;  % 45 degrees in radians
endAngle_true = 3*pi/4;  % 135 degrees in radians

theta_true = linspace(startAngle_true, endAngle_true, 100);
x_true = center_xCoor_true + radius_true * cos(theta_true);
x_true = x_true';
y_true = center_yCoor_true + radius_true * sin(theta_true);
y_true = y_true';

% Add some noise to the data
noise = 0.1 * randn(size(x_true,1),1);
x_noisy = x_true + noise;
y_noisy = y_true + noise;

% Initial parameter estimates (these are just initial guesses)
initial_center_xCoor_guess = 5;
initial_center_yCoor_guess = 5;
initial_radius_guess = 5;
initial_startAngle_guess = 1;
initial_endAngle_guess = 1;

% Combine noisy data into a single vector of points
noisyData = [x_noisy, y_noisy];

initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

[fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, 123);
disp(fittedArc);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Arc', 'Location', 'best');

%% Arc - Case 4

% Generate synthetic noisy data resembling a circular arc
center_xCoor_true = 20;
center_yCoor_true = 10;
radius_true = 10;
startAngle_true = pi;  % 45 degrees in radians
endAngle_true = -pi;  % 135 degrees in radians

theta_true = linspace(startAngle_true, endAngle_true, 100);
x_true = center_xCoor_true + radius_true * cos(theta_true);
x_true = x_true';
y_true = center_yCoor_true + radius_true * sin(theta_true);
y_true = y_true';

% Add some noise to the data
noise = 0.2 * randn(size(x_true,1),1);
x_noisy = x_true + noise;
y_noisy = y_true + noise;

% Initial parameter estimates (these are just initial guesses)
initial_center_xCoor_guess = 5;
initial_center_yCoor_guess = 5;
initial_radius_guess = 5;
initial_startAngle_guess = 1;
initial_endAngle_guess = 1;

% Combine noisy data into a single vector of points
noisyData = [x_noisy, y_noisy];

initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

[fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, 1234);
disp(fittedArc);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Arc', 'Location', 'best');