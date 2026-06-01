%% script_test_fcn_geometry_maxRadiusInsideBox
% Exercises the function: fcn_geometry_maxRadiusInsideBox

% 2024_04_14 - S. Brennan
% - wrote the code

close all;

%% Basic example - normal case
figNum = 1;
figure(figNum);
clf;

box_width  = 4;
box_height = 1;
max_radius = fcn_geometry_maxRadiusInsideBox(box_width, box_height, figNum);
assert(isequal(round(max_radius,4),round(2.5,4)));


%% Basic example - height larger than L/2
figNum = 2;
figure(figNum);
clf;

box_width  = 4;
box_height = 3;
max_radius = fcn_geometry_maxRadiusInsideBox(box_width, box_height, figNum);
assert(isequal(round(max_radius,4),round(2,4)));








