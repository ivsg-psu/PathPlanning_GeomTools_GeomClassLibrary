% script_test_fcn_geometry_fitSpiral.m
% tests fcn_geometry_fitSpiral.m

% Revision history
% 2023_10_19 - Aneesh Batchu
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

%% Archimedean Spiral - Case 1

% The spiral is used, according to ASAM OpenDRIVE, to describe the
% transition from a line element to an arc element

% to ensure results are repeatable within a MATLAB session
rng('default')

% Define the coordinates of the two points
FirstPoint = [1 2]; % End Point of a line
LastPoint = [3 4];  % Starting point of an arc

% Calculate the angles of the points
thetaFirstPoint = atan2(FirstPoint(1,2), FirstPoint(1,1));
thetaLastPoint = atan2(LastPoint(1,2), LastPoint(1,1));

% Calculate the values of a (Initial Radius) and b (Slope)
r_FirstPoint = sqrt(sum(FirstPoint.^2));
r_LastPoint = sqrt(sum(LastPoint.^2));
b = (r_LastPoint - r_FirstPoint) / (thetaLastPoint - thetaFirstPoint);
a = r_FirstPoint - b * thetaFirstPoint;

% Generate points along the fitted Archimedean spiral
theta = linspace(thetaFirstPoint, thetaLastPoint,50);
r = a + b * theta;
x = r .* cos(theta);
y = r .* sin(theta);
NoNoiseData = [x', y'];
% plot(x, y, 'blue', 'LineWidth', 2);
% hold on
% scatter(x, y,'red','filled','SizeData', 30)

% Adding noise to the points in between the end point of a line and the
% start point of an arc
data = [x(2:end-1)' y(2:end-1)'];
noise = 0.05*randn(size(data)); % noise

noisyData = data + noise;
noisyData = [x(1), y(1); noisyData; x(end), y(end)];

% initial guess for a (Initial radius) and b (Slope)
initialGuess = [1 1];
fig_num = 234;

[fittedSpiral,fittedParameters] = fcn_geometry_fitSpiral(noisyData, initialGuess, fig_num);

disp(fittedParameters)
disp(fittedSpiral)

% RMS Error of fitted Spiral data with the original data without noise
RMSE = sqrt(mean(sum((NoNoiseData-fittedSpiral).^2,2)));
disp('RMSE: ')
disp(RMSE)
 


%% Archimedean Spiral - Case 2

% Define the coordinates of the two points
FirstPoint = [1 1]; % End Point of a line
LastPoint = [6 4];  % Starting point of an arc

% Calculate the angles of the points
thetaFirstPoint = atan2(FirstPoint(1,2), FirstPoint(1,1));
thetaLastPoint = atan2(LastPoint(1,2), LastPoint(1,1));

% Calculate the values of a and b
r_FirstPoint = sqrt(sum(FirstPoint.^2));
r_LastPoint = sqrt(sum(LastPoint.^2));
b = (r_LastPoint - r_FirstPoint) / (thetaLastPoint - thetaFirstPoint);
a = r_FirstPoint - b * thetaFirstPoint;

% Generate points along the fitted Archimedean spiral
theta = linspace(thetaFirstPoint, thetaLastPoint, 50);
r = a + b * theta;
x = r .* cos(theta);
y = r .* sin(theta);

data = [x(2:end-1)' y(2:end-1)'];
noise = 0.02*randn(size(data));

noisyData = data + noise;
noisyData = [x(1), y(1); noisyData; x(end), y(end)];

initialGuess = [1 1];

[fittedSpiral,fittedParameters] = fcn_geometry_fitSpiral(noisyData, initialGuess, 12);

disp(fittedParameters)
disp(fittedSpiral)