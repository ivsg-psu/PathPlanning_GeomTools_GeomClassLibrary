% script_test_fcn_geometry_fitStraightLine.m
% tests fcn_geometry_fitStraightLine.m

% Revision history
% 2023_10_17 - Aneesh Batchu
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

%% Staright Line - Case 1

% Generate random x-values 
x = 15 * rand(100, 1);

% Create corresponding y-values with some noise
noise = 2 * randn(100, 1); 
slope = 8; 
intercept = 5; % Y-intercept
y = slope * x + intercept + noise;

noisyData = [x, y];

[fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData,1);
disp(fittedLine);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Line', 'Location', 'best');

%% Staright Line - Case 2

% Generate random x-values 
x = 15 * rand(100, 1);

% Create corresponding y-values with some noise
noise = 0.3* randn(100, 1); 
slope = 0; 
intercept = 5; % Y-intercept
y = slope * x + intercept + noise;

noisyData = [x, y];

[fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData,12);
disp(fittedLine);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Line', 'Location', 'best');

%% Staright Line - Case 3

% Generate random x-values 
x = 2 * rand(100, 1);

% Create corresponding y-values with some noise
noise = 0.1 * randn(100, 1); 
slope = 3; 
intercept = 0; % Y-intercept
y = slope * x + intercept + noise;

noisyData = [x, y];

[fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData,123);
disp(fittedLine);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Line', 'Location', 'best');

%% Staright Line - Case 3

% Generate random x-values 
x = 2 * rand(100, 1);

% Create corresponding y-values with some noise
noise = 0.1 * randn(100, 1); 
slope = 1; 
intercept = 0; % Y-intercept
y = slope * x + intercept + noise;

noisyData = [x, y];

[fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData,1234);
disp(fittedLine);
disp(fittedParameters);
legend('Noisy Points', 'Fitted Line', 'Location', 'best');
