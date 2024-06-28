%% script_test_fcn_geometry_fitRadiusToPoints
% Exercises the function: fcn_geometry_fitRadiusToPoints

% 2024_06_27 - S. Brennan
% -- wrote the code

close all;

%% Basic example - normal case
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
% arc_pattern = [...
%     1/20, 15; 
%     0 20;
%     -1/5 10; 
%     0 10;
%     1/15 40; 
%     0 15
%     -1/10 20];

arc_pattern = [...
    1/20, 45; 
];

M = 10; % How many points per meter
sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, fig_num);

% Call the function
[radius, radius_maximum] = fcn_geometry_fitRadiusToPoints(test_points, fig_num);

% Check sizes
assert(isequal(size(radius),[1 1]));
assert(isequal(size(radius_maximum),[1 1]));

% Check values
% assert(isequal(size(radius),[1 1]));




