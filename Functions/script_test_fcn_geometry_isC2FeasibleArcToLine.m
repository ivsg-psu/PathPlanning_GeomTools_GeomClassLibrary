%% script_test_fcn_geometry_isC2FeasibleArcToLine
% Tests the function fcn_geometry_isC2FeasibleArcToLine

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

%% 1.1 Basic test - feasible

fig_num = 1101;
figure(fig_num);
clf;

% Fill in line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_xy = [0 -1];
line_angle   = 0;

line_parameters(1,1:2) = line_base_xy;
line_parameters(1,3)   = line_angle;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_radius               = 1;
x_offset = 0; 
y_offset = 0; 
circle_center_xy            = [x_offset circle_radius+y_offset];

circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

% Set threshold and margin
threshold = 0;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_line_parameters] = fcn_geometry_isC2FeasibleArcToLine( circle_parameters, line_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_line_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-1));
assert(isequal(round(closest_feasible_line_parameters,4),[0    -1    0]));


%% 1.2 Basic test - NOT feasible

fig_num = 1201;
figure(fig_num);
clf;

% Fill in line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_xy = [0 0.4];
line_angle   = 0;

line_parameters(1,1:2) = line_base_xy;
line_parameters(1,3)   = line_angle;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_radius               = 1;
x_offset = 0; 
y_offset = 0; 
circle_center_xy            = [x_offset circle_radius+y_offset];

circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

% Set threshold and margin
threshold = 0.2;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_line_parameters] = fcn_geometry_isC2FeasibleArcToLine( circle_parameters, line_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_line_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.4));
assert(isequal(round(closest_feasible_line_parameters,4),[0    0    0]));

%%%% 
% Show that alignment is not possible with original parameters
circle_line_parameters = [line_base_xy inf];
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle_parameters, circle_line_parameters, [], []);


% Check results - show this is not feasible - either NaN or length is zero
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1,4),4),0));


%% 1.3 Basic test - marginally infeasible

fig_num = 1301;
figure(fig_num);
clf;

% Fill in line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_xy = [0 0.4];
line_angle   = 0;

line_parameters(1,1:2) = line_base_xy;
line_parameters(1,3)   = line_angle;

% Fill in circle2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_radius               = 1;
x_offset = 0; 
y_offset = 0; 
circle_center_xy            = [x_offset circle_radius+y_offset];

circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

% Set threshold and margin
threshold = 0.4;
in_boundary_margin = [];

% Call function
[flag_is_feasible, feasibility_distance, closest_feasible_line_parameters] = fcn_geometry_isC2FeasibleArcToLine( circle_parameters, line_parameters, (threshold), (in_boundary_margin), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));
assert(isequal(size(closest_feasible_line_parameters),[1 3]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),0.4));
assert(isequal(round(closest_feasible_line_parameters,4),[0    0     0]));

%%%% 
% Show that alignment is not possible with original parameters
circle_line_parameters = [line_base_xy inf];
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle_parameters, circle_line_parameters, [], []);


% Check results - show this is not feasible - either NaN or length is zero
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1,4),4),0));



%%%% 
% Show that alignment again is not possible using modified parameters that
% are exactly on feasible boundary

circle_line_parameters = [closest_feasible_line_parameters(1,1:2) inf];
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle_parameters, circle_line_parameters, [], []);


% Check results - show this is not feasible - either NaN or length is zero
assert(all(isnan(spiral_join_parameters)) || isequal(round(spiral_join_parameters(1,4),4),0));


%%%
% Call function with a user-defined margin that is small

% Set margin to push solution deeper into feasible area
in_boundary_margin = 0.01; % Units are meters
[flag_is_feasible, feasibility_distance, closest_feasible_line_parameters] = fcn_geometry_isC2FeasibleArcToLine( circle_parameters, line_parameters, (threshold), (in_boundary_margin), (fig_num));

assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),0.4));
assert(isequal(round(closest_feasible_line_parameters,4),[0    -0.01    0]));

%%%% 
% Show that alignment again now possible with modified parameters that
% are inside feasible boundary

circle_line_parameters = [closest_feasible_line_parameters(1,1:2) inf];
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle_parameters, circle_line_parameters, [], 234456);

% Check results - show this is feasible
assert(~all(isnan(spiral_join_parameters)));

