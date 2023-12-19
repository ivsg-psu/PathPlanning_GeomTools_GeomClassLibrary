% script_test_fcn_geometry_circleHoughFit
% Exercises the function: fcn_geometry_circleHoughFit

% Revision history:
% 2023_12_15
% -- wrote the code

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Fill test data - 1 segment
fig_num = 23;
% figure(fig_num);
% clf;
% hold on;
% axis equal
% grid on;

seed_points = [2 3; 4 5; 6 3]; %; 1 1];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

%% Example - 1
fig_num = 1;

inputPoints = test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.1;

[fittedParameters, agreementIndices] = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, fig_num);

agreements = sum(agreementIndices,2);
% max(agreements)
