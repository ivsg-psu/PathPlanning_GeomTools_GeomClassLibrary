%% script_test_fcn_geometry_gapCircleToCircle
% Tests the function: fcn_geometry_gapCircleToCircle

% Revision history:
% 2024_05_06 - S. Brennan
% - wrote the code

close all

%% Basic test - circle 1 and circle 2 are large and small, counter-clockwise
% If larger to smaller, then small one has to be inside. If any outside, it
% will not join.
figNum = 1;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 1.1];
flag_circle2_is_inside = 1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),0.0895));

%% Basic test - circle 1 and circle 2 are large and small, counter-clockwise
% show that if circle is outside, it will not join
figNum = 2;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 0.8];
flag_circle2_is_inside = 1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),-0.2091));



%% Basic test - circle 1 and circle 2 are small and large, counter-clockwise
% If smaller to larger, then small one has to be inside. If any outside, it
% will not join.
figNum = 3;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 1.8];
flag_circle2_is_inside = 1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),0.1754));

%% Basic test - circle 1 and circle 2 are small and large, counter-clockwise
% Show that if circle is outside, it will not join
figNum = 4;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 2.2];
flag_circle2_is_inside = 1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),-0.2166));



%% Basic test - circle 1 and circle 2 are large and small, CLOCKWISE
% If clockwise, the small circle must be completely outside the large
% circle
figNum = 5;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 -1.1];
flag_circle2_is_inside = -1;


gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),0.1049));

%% Basic test - circle 1 and circle 2 are large and small, CLOCKWISE
% show that if circle is outside, it will not join
figNum = 6;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 -0.9];
flag_circle2_is_inside = -1;


gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),-0.0949));



%% Basic test - circle 1 and circle 2 are small and large, CLOCKWISE
% If smaller to larger, then small one has to be inside. If any outside, it
% will not join.
figNum = 7;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 -2.2];
flag_circle2_is_inside = -1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);

% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),0.2062));

%% Basic test - circle 1 and circle 2 are small and large, CLOCKWISE
% Show that if circle is outside, it will not join
figNum = 8;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 -1.8];
flag_circle2_is_inside = -1;


gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);


% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),-0.1929));


%% FAIL test - circle 1 and circle 2 are small and large, COUNTERCLOCKWISE
% Show that if circle is outside, it will not join
figNum = 9999;
figure(figNum); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 3;
circle2_center_XY = [0 2.999];
flag_circle2_is_inside = 1;

gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_XY, flag_circle2_is_inside, figNum);



% Check size of results
assert(isequal(size(gap_between_circles),[1 1]));

% Check results
assert(isequal(round(gap_between_circles,4),0.0010));
