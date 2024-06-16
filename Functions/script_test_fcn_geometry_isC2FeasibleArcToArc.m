%% script_test_fcn_geometry_isC2FeasibleArcToArc
% Tests the function fcn_geometry_isC2FeasibleArcToArc

% Revision history:
% 2024_05_26 - S Brennan
% -- wrote the code

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
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

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
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

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
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

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
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

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
circle1_radius               = 2;
x_offset = 0; 
y_offset = 0; 
circle1_center_xy            = [x_offset circle1_radius+y_offset];

circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_radius               = 3;
x_offset = 0; 
y_offset = -0.1; 
circle2_center_xy            = [x_offset -1*(circle2_radius+y_offset)];

circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_arc2_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[0   -2.9500    2.9500]));

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
[flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(circle1_parameters, circle2_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));
assert(isequal(round(closest_feasible_arc2_parameters,4),[ 0   -2.9571    2.9429]));

%%%% 
% Show that alignment again is possible with modified parameters that
% are inside feasible boundary
spiral_fig_num = 6547;
flag_circle2_is_counterclockwise = 0;
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, closest_feasible_arc2_parameters, flag_circle2_is_counterclockwise, spiral_fig_num);

% Check results - show this is not feasible
assert(~all(isnan(spiral_join_parameters)));





