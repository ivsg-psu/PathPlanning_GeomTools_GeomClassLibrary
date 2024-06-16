%% script_test_fcn_geometry_spiralFromCircleToCircle
% Tests the function: fcn_geometry_spiralFromCircleToCircle

% Revision history:
% 2024_04_24 - S. Brennan
% -- wrote the code
% 2024_05_28 - S. Brennan
% -- modified for parameter inputs, instead of raw inputs

close all

%% Basic test - circle 1 and circle 2 are large and small, counter-clockwise
% If larger to smaller, then small one has to be inside. If any outside, it
% will not join.
fig_num = 1;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 1.1];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[ -0.5746    0.0555   -0.1927    1.8195    0.3333    1.0000]));

%% Basic test - circle 1 and circle 2 are large and small, counter-clockwise
% show that if circle is outside, it will not join
fig_num = 2;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 0.8];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];


spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(all(isnan(spiral_join_parameters)));



%% Basic test - circle 1 and circle 2 are small and large, counter-clockwise
% If smaller to larger, then small one has to be inside. If any outside, it
% will not join.
fig_num = 3;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 1.8];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[-0.9770    1.2132   -1.7856    3.0183    1.0000    0.5000]));

%% Basic test - circle 1 and circle 2 are small and large, counter-clockwise
% Show that if circle is outside, it will not join
fig_num = 4;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 2.2];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(all(isnan(spiral_join_parameters)));



%% Basic test - circle 1 and circle 2 are large and small, CLOCKWISE
% If clockwise, the small circle must be completely outside the large
% circle
fig_num = 5;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 -1.1];
flag_circle2_is_counterclockwise = -1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];


spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[-0.5312    0.0474   -0.1780    1.3893    0.3333   -1.0000]));

%% Basic test - circle 1 and circle 2 are large and small, CLOCKWISE
% show that if circle is outside, it will not join
fig_num = 6;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 3;
circle2_radius = 1;
circle2_center_XY = [0.2 -0.9];
flag_circle2_is_counterclockwise = -1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(all(isnan(spiral_join_parameters)));



%% Basic test - circle 1 and circle 2 are small and large, CLOCKWISE
% If smaller to larger, then small one has to be inside. If any outside, it
% will not join.
fig_num = 7;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 -2.2];
flag_circle2_is_counterclockwise = -1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);

% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[-0.7503    0.3388   -0.8484    1.8584    1.0000   -0.5000]));

%% Basic test - circle 1 and circle 2 are small and large, CLOCKWISE
% Show that if circle is outside, it will not join
fig_num = 8;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 2;
circle2_center_XY = [0.2 -1.8];
flag_circle2_is_counterclockwise = -1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];


[spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);


% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(all(isnan(spiral_join_parameters)));
assert(space_between_circles<0);

%% FAIL test - circle 1 and circle 2 are small and large, COUNTERCLOCKWISE
% Show that if circle is outside, it will not join
fig_num = 9999;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = 3;
circle2_center_XY = [0 2.999];
flag_circle2_is_counterclockwise = 1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

[spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);


% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[-0.0947    0.0045   -0.0949    0.1897    1.0000    0.3333]));
assert(space_between_circles>0);

%% FAIL test - not sure why this is not working well?
% Show that if circle is outside, it will not join
fig_num = 99988;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 162.3452;
circle2_radius = 367.6955;
circle2_center_XY = [-36.3698  364.4481];
flag_circle2_is_counterclockwise = 1;

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

[spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);


% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[ 27.4677    2.3405    0.1700    2.6105    0.0062    0.0027]));
assert(space_between_circles>0);

%% Circle to line tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _          _             _          _      _
%  / ____(_)        | |           | |        | |    (_)
% | |     _ _ __ ___| | ___  ___  | |_ ___   | |     _ _ __   ___  ___
% | |    | | '__/ __| |/ _ \/ __| | __/ _ \  | |    | | '_ \ / _ \/ __|
% | |____| | | | (__| |  __/\__ \ | || (_) | | |____| | | | |  __/\__ \
%  \_____|_|_|  \___|_|\___||___/  \__\___/  |______|_|_| |_|\___||___/
%
%
% http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circles%20to%20Lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Line test - circle 1 and line, line below x-axis - feasible
% Show that if line is outside, it will join
fig_num = 101;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = inf;
circle2_center_XY = [0 -0.2];
flag_circle2_is_counterclockwise = [];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

[spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);


% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(isequal(round(spiral_join_parameters,4),[-0.9002    0.5645   -1.1202    2.2403    1.0000         0]));
assert(space_between_circles==0.2);


%% Line test - circle 1 and line, line above x-axis - not feasible
% Show that if line is not outside, it will not join
fig_num = 102;
figure(fig_num); clf;


% Try positive curvature
circle1_radius = 1;
circle2_radius = inf;
circle2_center_XY = [0 0.2];
flag_circle2_is_counterclockwise = [];

circle1_parameters = [0 circle1_radius  circle1_radius];
circle2_parameters = [circle2_center_XY circle2_radius];

[spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, flag_circle2_is_counterclockwise, fig_num);


% Check size of results
assert(isequal(size(spiral_join_parameters),[1 6]));

% Check results
assert(all(isnan(spiral_join_parameters)));
assert(space_between_circles<0);