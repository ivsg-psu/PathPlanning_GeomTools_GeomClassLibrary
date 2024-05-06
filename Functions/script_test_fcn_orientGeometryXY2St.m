%% script_test_fcn_orientGeometryXY2St
% Exercises the function: fcn_orientGeometryXY2St

% Revision history:
% 2024_05_02 - S. Brennan
% -- wrote the code

close all;

%% Basic test 1.1 - a line segment alone
fig_num = 11;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_angle = 30*pi/180;
segment_base_point_xy = [ 2 3];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 3;
segment_s_end   = 7;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

primary_parameters_type_string = 'segment';
primary_parameters = segment_parameters;
secondary_parameters_type_strings = [];
secondary_parameters = [];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 6]));
assert(isequal(size(st_secondary_parameters),[1 1]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[1.0000    0.0000   -4.0000   -0.0000         0    4.0000]));
assert(isempty(st_secondary_parameters{1}));
assert(dist(St_transform,se2([0.8660    0.5000  -10.2321;    -0.5000    0.8660   -1.5981;   0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle*180/pi,4),-30));
assert(isequal(flag_primary_parameter_is_flipped,0));

%% Basic test 1.21 - an arc alone, counter-clockwise
fig_num = 121;
figure(fig_num); clf;

arc_center_xy            = [2 3];
arc_radius               = 4;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

primary_parameters_type_string = 'arc';
primary_parameters = arc_parameters;
secondary_parameters_type_strings = [];
secondary_parameters = [];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 7]));
assert(isequal(size(st_secondary_parameters),[1 1]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[0    4.0000    4.0000   -2.7925   -1.5708         0    1.0000]));
assert(isempty(st_secondary_parameters{1}));
assert(dist(St_transform,se2([    0.0000    1.0000   -3.0000;    -1.0000    0.0000    6.0000;  0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle,4),-1.5708));
assert(isequal(flag_primary_parameter_is_flipped,0));


%% Basic test 1.22 - an arc alone, clockwise
fig_num = 122;
figure(fig_num); clf;

arc_center_xy            = [2 3];
arc_radius               = 4;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

primary_parameters_type_string = 'arc';
primary_parameters = arc_parameters;
secondary_parameters_type_strings = [];
secondary_parameters = [];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 7]));
assert(isequal(size(st_secondary_parameters),[1 1]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[0    4.0000    4.0000   -0.3491   -1.5708         0    1.0000]));
assert(isempty(st_secondary_parameters{1}));
assert(dist(St_transform,se2([-0.0000   -1.0000   3.0000;    -1.0000    0.0000   6.0000;  0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle,4),1.5708));
assert(isequal(flag_primary_parameter_is_flipped,1));


%% Basic test 2.1 - a line segment with other geometries
fig_num = 21;
figure(fig_num); clf;
clear secondary_parameters secondary_parameters_type_strings

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_angle = 30*pi/180;
segment_base_point_xy = [ 2 3];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 3;
segment_s_end   = 7;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

primary_parameters_type_string = 'segment';
primary_parameters = segment_parameters;


%%%%
% Fill in the other geometries - one of each type

% A test line
segment_angle = 60*pi/180;
segment_base_point_xy = [ 6 7];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 1;
segment_s_end   = 4;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

secondary_parameters_type_strings{1} = 'segment';
secondary_parameters{1}              = segment_parameters;

% A test arc counter-clockwise
arc_center_xy            = [1 2];
arc_radius               = 2;
arc_vector_start         = [cos(30*pi/180) sin(30*pi/180)];
arc_vector_end           = [cos(80*pi/180) sin(80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{2} = 'arc';
secondary_parameters{2}              = arc_parameters;

% A test arc clockwise
arc_center_xy            = [4 7];
arc_radius               = 1;
arc_vector_start         = [cos(-30*pi/180) sin(-30*pi/180)];
arc_vector_end           = [cos(-80*pi/180) sin(-80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{3} = 'arc';
secondary_parameters{3}              = arc_parameters;

% A test line
% [unit_projection_vector_x,
%     unit_projection_vector_y,
%     base_point_x,
%     base_point_y,
%     ]
secondary_parameters_type_strings{4} = 'line';
secondary_parameters{4}              = [cos(-30*pi/180) sin(-30*pi/180) 2 2];

% A test circle
% [circleCenter_x.
%     circleCenter_y,
%     radius]
secondary_parameters_type_strings{5} = 'circle';
secondary_parameters{5}              = [-1 5 2];

% A test spiral
%  [spiralLength,  % the s-coordinate length allowed
%   h0,  % The initial heading
%   x0,  % The initial x value
%   y0,  % The initial y value
%   K0,  % The initial curvature
%   Kf   % The final curvature
% ]
secondary_parameters_type_strings{6} = 'spiral';
secondary_parameters{6}              = [4 -40*pi/180 5 5 -1 4];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 6]));
assert(isequal(size(st_secondary_parameters),[1 6]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[1.0000    0.0000   -4.0000   -0.0000         0    4.0000]));
assert(isequal(round(st_secondary_parameters{1},4), [0.8660    0.5000   -0.6699    1.9641         0    3.0000]));
assert(isequal(round(st_secondary_parameters{2},4), [-8.3660   -0.3660    2.0000         0    0.8727         0    1.0000]));
assert(isequal(round(st_secondary_parameters{3},4), [-3.2679    2.4641    1.0000   -1.0472   -1.9199         0         0]));
assert(isequal(round(st_secondary_parameters{4},4), [0.5000   -0.8660   -7.5000   -0.8660]));
assert(isequal(round(st_secondary_parameters{5},4), [-8.5981    3.2321    2.0000]));
assert(isequal(round(st_secondary_parameters{6},4), [4.0000   -1.2217   -3.4019    0.2321   -1.0000    4.0000]));

assert(dist(St_transform,se2([0.8660    0.5000  -10.2321;    -0.5000    0.8660   -1.5981;   0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle*180/pi,4),-30));
assert(isequal(flag_primary_parameter_is_flipped,0));

%% Basic test 2.21 - an arc, counter-clockwise with other geometries
fig_num = 221;
figure(fig_num); clf;
clear secondary_parameters secondary_parameters_type_strings

% Fill in the arc parameters
arc_center_xy            = [2 3];
arc_radius               = 4;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

primary_parameters_type_string = 'arc';
primary_parameters = arc_parameters;


%%%%
% Fill in the other geometries - one of each type

% A test line
segment_angle = 60*pi/180;
segment_base_point_xy = [ 6 7];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 1;
segment_s_end   = 4;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

secondary_parameters_type_strings{1} = 'segment';
secondary_parameters{1}              = segment_parameters;

% A test arc counter-clockwise
arc_center_xy            = [1 2];
arc_radius               = 2;
arc_vector_start         = [cos(30*pi/180) sin(30*pi/180)];
arc_vector_end           = [cos(80*pi/180) sin(80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{2} = 'arc';
secondary_parameters{2}              = arc_parameters;

% A test arc clockwise
arc_center_xy            = [4 7];
arc_radius               = 1;
arc_vector_start         = [cos(-30*pi/180) sin(-30*pi/180)];
arc_vector_end           = [cos(-80*pi/180) sin(-80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{3} = 'arc';
secondary_parameters{3}              = arc_parameters;

% A test line
% [unit_projection_vector_x,
%     unit_projection_vector_y,
%     base_point_x,
%     base_point_y,
%     ]
secondary_parameters_type_strings{4} = 'line';
secondary_parameters{4}              = [cos(-30*pi/180) sin(-30*pi/180) 2 2];

% A test circle
% [circleCenter_x.
%     circleCenter_y,
%     radius]
secondary_parameters_type_strings{5} = 'circle';
secondary_parameters{5}              = [-1 5 2];

% A test spiral
%  [spiralLength,  % the s-coordinate length allowed
%   h0,  % The initial heading
%   x0,  % The initial x value
%   y0,  % The initial y value
%   K0,  % The initial curvature
%   Kf   % The final curvature
% ]
secondary_parameters_type_strings{6} = 'spiral';
secondary_parameters{6}              = [4 -40*pi/180 5 5 -1 4];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 7]));
assert(isequal(size(st_secondary_parameters),[1 6]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[0    4.0000    4.0000   -2.7925   -1.5708         0    1.0000]));
assert(isequal(round(st_secondary_parameters{1},4), [0.8660   -0.5000    4.8660   -0.5000         0    3.0000]));
assert(isequal(round(st_secondary_parameters{2},4), [-1.0000    5.0000    2.0000   -1.0472   -0.1745         0    1.0000]));
assert(isequal(round(st_secondary_parameters{3},4), [4.0000    2.0000    1.0000   -2.0944   -2.9671         0         0]));
assert(isequal(round(st_secondary_parameters{4},4), [-0.5000   -0.8660   -1.0000    4.0000]));
assert(isequal(round(st_secondary_parameters{5},4), [2 7 2]));
assert(isequal(round(st_secondary_parameters{6},4), [4.0000   -2.2689    2.0000    1.0000   -1.0000    4.0000]));

assert(dist(St_transform,se2([0.0000    1.0000   -3.0000;    -1.0000    0.0000    6.0000;   0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle*180/pi,4),-90));
assert(isequal(flag_primary_parameter_is_flipped,0));

%% Basic test 2.22 - an arc,clockwise with other geometries
fig_num = 222;
figure(fig_num); clf;
clear secondary_parameters secondary_parameters_type_strings

% Fill in the arc parameters
arc_center_xy            = [2 3];
arc_radius               = 4;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

primary_parameters_type_string = 'arc';
primary_parameters = arc_parameters;


%%%%
% Fill in the other geometries - one of each type

% A test line
segment_angle = 60*pi/180;
segment_base_point_xy = [ 6 7];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 1;
segment_s_end   = 4;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

secondary_parameters_type_strings{1} = 'segment';
secondary_parameters{1}              = segment_parameters;

% A test arc counter-clockwise
arc_center_xy            = [1 2];
arc_radius               = 2;
arc_vector_start         = [cos(30*pi/180) sin(30*pi/180)];
arc_vector_end           = [cos(80*pi/180) sin(80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{2} = 'arc';
secondary_parameters{2}              = arc_parameters;

% A test arc clockwise
arc_center_xy            = [4 7];
arc_radius               = 1;
arc_vector_start         = [cos(-30*pi/180) sin(-30*pi/180)];
arc_vector_end           = [cos(-80*pi/180) sin(-80*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;

secondary_parameters_type_strings{3} = 'arc';
secondary_parameters{3}              = arc_parameters;

% A test line
% [unit_projection_vector_x,
%     unit_projection_vector_y,
%     base_point_x,
%     base_point_y,
%     ]
secondary_parameters_type_strings{4} = 'line';
secondary_parameters{4}              = [cos(-30*pi/180) sin(-30*pi/180) 2 2];

% A test circle
% [circleCenter_x.
%     circleCenter_y,
%     radius]
secondary_parameters_type_strings{5} = 'circle';
secondary_parameters{5}              = [-1 5 2];

% A test spiral
%  [spiralLength,  % the s-coordinate length allowed
%   h0,  % The initial heading
%   x0,  % The initial x value
%   y0,  % The initial y value
%   K0,  % The initial curvature
%   Kf   % The final curvature
% ]
secondary_parameters_type_strings{6} = 'spiral';
secondary_parameters{6}              = [4 -40*pi/180 5 5 -1 4];

% Call the function
[st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num));

% Check size of results
assert(isequal(size(st_primary_parameters),[1 7]));
assert(isequal(size(st_secondary_parameters),[1 6]));
assert(isequal(size(St_transform),[1 1]));
assert(isequal(size(rotation_angle),[1 1]));
assert(isequal(size(flag_primary_parameter_is_flipped),[1 1]));


% Check results
assert(isequal(round(st_primary_parameters,4),[0    4.0000    4.0000   -0.3491   -1.5708         0    1.0000]));
assert(isequal(round(st_secondary_parameters{1},4), [-0.8660   -0.5000   -4.8660   -0.5000         0    3.0000]));
assert(isequal(round(st_secondary_parameters{2},4), [1.0000    5.0000    2.0000   -2.0944   -2.9671         0         0]));
assert(isequal(round(st_secondary_parameters{3},4), [-4.0000    2.0000    1.0000   -1.0472   -0.1745         0    1.0000]));
assert(isequal(round(st_secondary_parameters{4},4), [0.5000   -0.8660    1.0000    4.0000]));
assert(isequal(round(st_secondary_parameters{5},4), [-2 7 2]));
assert(isequal(round(st_secondary_parameters{6},4), [4.0000   -0.8727   -2.0000    1.0000    1.0000   -4.0000]));

assert(dist(St_transform,se2([0.0000    -1.0000   3.0000;    -1.0000    0.0000    6.0000;   0         0    1.0000]))<0.001);
assert(isequal(round(rotation_angle*180/pi,4),90));
assert(isequal(flag_primary_parameter_is_flipped,1));

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_orientGeometryXY2St(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end