%% script_test_fcn_geometry_curvatureAtPoint
% Exercises the function: fcn_geometry_curvatureAtPoint

% 2024_06_27 - S. Brennan
% - wrote the code

% close all;

%% Basic example - normal case in arc
figNum = 1;
figure(figNum);
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
    0 10;
    1/20, 45; 
    0 20;
];

M = 10; % How many points per meter
sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

[points_to_fit, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

index_to_test = 200;
[point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, index_to_test, [], (figNum));

% Check sizes
assert(isequal(size(point_curvature),[1 1]));
assert(isequal(size(point_circle_center),[1 2]));
assert(isequal(size(index_range),[1 1]));
assert(isequal(size(point_curvature_minimum),[1 1]));


%% Basic example - normal case on line portion
figNum = 2;
figure(figNum);
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
    0 10;
    1/20, 45; 
    0 20;
];

M = 10; % How many points per meter
sigma = 0.05; % The standard deviation in the points relative to the perfect function fit, in meters

[points_to_fit, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

index_to_test = 600;
[point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, index_to_test, [], (figNum));

% Check sizes
assert(isequal(size(point_curvature),[1 1]));
assert(isequal(size(point_circle_center),[1 2]));
assert(isequal(size(index_range),[1 1]));
assert(isequal(size(point_curvature_minimum),[1 1]));

%% Basic example - normal case near end
figNum = 3;
figure(figNum);
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
    0 10;
    1/20, 45; 
    0 20;
];

M = 10; % How many points per meter
sigma = 0.05; % The standard deviation in the points relative to the perfect function fit, in meters

[points_to_fit, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

index_to_test = 20;
[point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, index_to_test, [], (figNum));

% Check sizes
assert(isequal(size(point_curvature),[1 1]));
assert(isequal(size(point_circle_center),[1 2]));
assert(isequal(size(index_range),[1 1]));
assert(isequal(size(point_curvature_minimum),[1 1]));

