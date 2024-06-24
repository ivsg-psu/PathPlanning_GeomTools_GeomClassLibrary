%% script_test_fcn_geometry_isC2FeasibleLineToArc
% Tests the function fcn_geometry_isC2FeasibleLineToArc

% Revision history:
% 2024_06_24 - S Brennan
% -- wrote the code

close all

% Sections:
% 1 - illustrative test cases
% 2 - feasible with 0 tolerance across each of 3 types
% 3 - infeasible with 0 tolerance across each of 3 types
% 4 - feasible with non-zero tolerance across each of 3 types
% 5 - infeasible with non-zero tolerance across each of 3 types

%% Illustrative test cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _____ _ _           _             _   _             _______        _      _____
% |_   _| | |         | |           | | (_)           |__   __|      | |    / ____|
%   | | | | |_   _ ___| |_ _ __ __ _| |_ ___   _____     | | ___  ___| |_  | |     __ _ ___  ___  ___
%   | | | | | | | / __| __| '__/ _` | __| \ \ / / _ \    | |/ _ \/ __| __| | |    / _` / __|/ _ \/ __|
%  _| |_| | | |_| \__ \ |_| | | (_| | |_| |\ V /  __/    | |  __/\__ \ |_  | |___| (_| \__ \  __/\__ \
% |_____|_|_|\__,_|___/\__|_|  \__,_|\__|_| \_/ \___|    |_|\___||___/\__|  \_____\__,_|___/\___||___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Illustrative%20Test%20Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1 Basic test - circle2 in circle1

fig_num = 1101;
figure(fig_num);
clf;

% Fill in circle1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
circle1_center_xy            = [x_offset circle1_radius+y_offset];

circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_radius               = 0.5;
x_offset = 0; 
y_offset = -0.1; 
circle2_center_xy            = [x_offset circle2_radius+y_offset];

circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    0.4500    0.4500]));

%%%% 
% Show that alignment is not possible with original
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    0.4571    0.4429]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 5858;
flag_circle2_is_counterclockwise = 1;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, flag_circle2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));



%% 1.2 Basic test - circle1 in circle2

fig_num = 1201;
figure(fig_num);
clf;

% Fill in circle1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_radius               = 2;
x_offset = 0; 
y_offset = 0; 
circle1_center_xy            = [x_offset circle1_radius+y_offset];

circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_radius               = 3;
x_offset = 0; 
y_offset = 0.1; 
circle2_center_xy            = [x_offset circle2_radius+y_offset];

circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    3.0500    3.0500]));

%%%% 
% Show that alignment is NOT possible with original parameters
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible. Either NaN values produced, or
% zero spiral length (the first parameter in spiral_join_parameters)
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1),4),0));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0    3.0429    3.0571]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 5478;
flag_circle2_is_counterclockwise = 1;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, flag_circle2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));



%% 1.3 Basic test - circle1 outside of circle2

fig_num = 1301;
figure(fig_num);
clf;

% Fill in circle1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_radius               = 4;
x_offset = 0; 
y_offset = 0; 
circle1_center_xy            = [x_offset circle1_radius+y_offset];

circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_radius               = 2;
x_offset = 0; 
y_offset = -0.1; 
circle2_center_xy            = [x_offset -1*(circle2_radius+y_offset)];

circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0   -1.9500    1.9500]));

%%%% 
% Show that alignment is NOT possible with original parameters
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible. Either NaN values produced, or
% zero spiral length (the first parameter in spiral_join_parameters)
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1),4),0));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0   -1.9571    1.9429]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 6547;
flag_circle2_is_counterclockwise = 0;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, flag_circle2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));


%% 2.1 Basic test - arc2 in arc1

fig_num = 2101;
figure(fig_num);
clf;

% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc1_radius               = 1;
arc1_x_offset             = 0; 
arc1_y_offset             = 0; 
arc1_start_angle_degrees  = -180;
arc1_end_angle_degrees    = -90;
arc1_is_counter_clockwise = 1;

arc1_center_xy            = [arc1_x_offset arc1_radius+arc1_y_offset];
arc1_vector_start         = [cos(arc1_start_angle_degrees*pi/180) sin(arc1_start_angle_degrees*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle_degrees  *pi/180) sin(arc1_end_angle_degrees  *pi/180)];
arc1_is_circle            = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc2_radius               = 0.5;
arc2_x_offset             = 0; 
arc2_y_offset             = -0.1; 
arc2_start_angle_degrees  = -90;
arc2_end_angle_degrees    = 10;
arc2_is_counter_clockwise = 1;

arc2_center_xy            = [arc2_x_offset arc2_radius+arc2_y_offset];
arc2_vector_start         = [cos(arc2_start_angle_degrees*pi/180) sin(arc2_start_angle_degrees*pi/180)];
arc2_vector_end           = [cos(arc2_end_angle_degrees  *pi/180) sin(arc2_end_angle_degrees  *pi/180)];
arc2_is_circle            = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 7]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    0.4500    0.4500 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment is not possible with original
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, arc2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    0.4571    0.4429 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 5858;
flag_arc2_is_counterclockwise = 1;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, flag_arc2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));



%% 2.2 Basic test - arc1 in arc2

fig_num = 2201;
figure(fig_num);
clf;

% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc1_radius               = 2;
arc1_x_offset             = 0; 
arc1_y_offset             = 0; 
arc1_start_angle_degrees  = -180;
arc1_end_angle_degrees    = -90;
arc1_is_counter_clockwise = 1;

arc1_center_xy            = [arc1_x_offset arc1_radius+arc1_y_offset];
arc1_vector_start         = [cos(arc1_start_angle_degrees*pi/180) sin(arc1_start_angle_degrees*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle_degrees  *pi/180) sin(arc1_end_angle_degrees  *pi/180)];
arc1_is_circle            = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc2_radius               = 3;
arc2_x_offset             = 0; 
arc2_y_offset             = 0.1; 
arc2_start_angle_degrees  = -90;
arc2_end_angle_degrees    = 10;
arc2_is_counter_clockwise = 1;

arc2_center_xy            = [arc2_x_offset arc2_radius+arc2_y_offset];
arc2_vector_start         = [cos(arc2_start_angle_degrees*pi/180) sin(arc2_start_angle_degrees*pi/180)];
arc2_vector_end           = [cos(arc2_end_angle_degrees  *pi/180) sin(arc2_end_angle_degrees  *pi/180)];
arc2_is_circle            = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 7]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0    3.0500    3.0500 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment is NOT possible with original parameters
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, arc2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible. Either NaN values produced, or
% zero spiral length (the first parameter in spiral_join_parameters)
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1),4),0));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0    3.0429    3.0571 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 5478;
flag_arc2_is_counterclockwise = 1;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, flag_arc2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));



%% 2.3 Basic test - arc1 outside of arc2

fig_num = 2301;
figure(fig_num);
clf;

% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc1_radius               = 4;
arc1_x_offset             = 0; 
arc1_y_offset             = 0; 
arc1_start_angle_degrees  = -180;
arc1_end_angle_degrees    = -90;
arc1_is_counter_clockwise = 1;

arc1_center_xy            = [arc1_x_offset arc1_radius+arc1_y_offset];
arc1_vector_start         = [cos(arc1_start_angle_degrees*pi/180) sin(arc1_start_angle_degrees*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle_degrees  *pi/180) sin(arc1_end_angle_degrees  *pi/180)];
arc1_is_circle            = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure

% Fill in arc 2
arc2_radius               = 2;
arc2_x_offset             = 0; 
arc2_y_offset             = -0.1; 
arc2_start_angle_degrees  = 110;
arc2_end_angle_degrees    = 10;
arc2_is_counter_clockwise = 0;

arc2_center_xy            = [arc2_x_offset -1*(arc2_radius+arc2_y_offset)];
arc2_vector_start         = [cos(arc2_start_angle_degrees*pi/180) sin(arc2_start_angle_degrees*pi/180)];
arc2_vector_end           = [cos(arc2_end_angle_degrees  *pi/180) sin(arc2_end_angle_degrees  *pi/180)];
arc2_is_circle            = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 7]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0   -1.9500    1.9500 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment is NOT possible with original parameters
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, arc2_parameters, [], []);

% Check results - show this is not feasible
assert(all(isnan(spiral_join_parameters)));



%%%% 
% Show that alignment again is not possible with modified parameters that
% are exactly on feasible boundary
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, [], []);

% Check results - show this is not feasible. Either NaN values produced, or
% zero spiral length (the first parameter in spiral_join_parameters)
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1),4),0));

%%%
% Call function with a user-defined margin that is small

% Set threshold and margin
threshold = 0;
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleLineToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0   -1.9571    1.9429 round(arc2_parameters(1,4:7),4)]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 6547;
flag_arc2_is_counterclockwise = 0;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(arc1_parameters, closest_feasible_arc2_parameters, flag_arc2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));




