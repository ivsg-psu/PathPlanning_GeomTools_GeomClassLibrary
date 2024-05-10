%% script_test_fcn_geometry_orientGeometrySt2XY
% Exercises the function: fcn_geometry_orientGeometrySt2XY

% Revision history:
% 2024_05_02 - S. Brennan
% -- wrote the code
% 2024_05_09 - S. Brennan
% -- fixed bug in segment calculation wherein unit vector gives NaN if
% start and end points are same

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

% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(isequal(size(XY_parameters{1}),[1 6]));

% Check results
assert(isequal(round(XY_parameters{1},4),[0.8660    0.5000    4.5981    4.5000         0    4.0000]));

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters{1};
% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, original_parameteters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results
assert(isequal(round(XY_parameters2{1},4),round(original_parameteters,4)));

%% Basic test 1.11 - a line segment of zero length
fig_num = 111;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_angle = 30*pi/180;
segment_base_point_xy = [ 2 3];
segment_unit_vector = [cos(segment_angle) sin(segment_angle)];
segment_s_start = 3;
segment_s_end   = 3;

segment_parameters(1,1:2) = segment_unit_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

primary_parameters_type_string = 'segment';
primary_parameters = segment_parameters;
secondary_parameters_type_strings = [];
secondary_parameters = [];

% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(isequal(size(XY_parameters{1}),[1 6]));

% Check results
assert(isequal(round(XY_parameters{1},4),[0.8660    0.5000    4.5981    4.5000         0    0]));

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters{1};
% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, original_parameteters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results
assert(isequal(round(XY_parameters2{1},4),round(original_parameteters,4)));

%% Basic test 1.12 - a NaN line segment
fig_num = 112;
figure(fig_num); clf;

segment_parameters   = nan(1,6);

primary_parameters_type_string = 'segment';
primary_parameters = segment_parameters;
secondary_parameters_type_strings = [];
secondary_parameters = [];

% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(isequal(size(XY_parameters{1}),[1 6]));

% Check results
assert(all(isnan(XY_parameters{1})));

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters{1};
% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, original_parameteters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results
assert(all(isnan(XY_parameters2{1})));


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


% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(isequal(size(XY_parameters{1}),[1 7]));

% Check results
assert(isequal(round(XY_parameters{1},4),[2.0000    3.0000    4.0000   -1.2217         0         0    1.0000]));

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters{1};
% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, original_parameteters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results
assert(isequal(round(XY_parameters2{1},4),round(original_parameteters,4)));



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

% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(isequal(size(XY_parameters{1}),[1 7]));

% Check results
assert(isequal(round(XY_parameters{1},4),[2.0000    3.0000    4.0000   -1.2217         0         0         0]));

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters{1};
% Call the function to convert from XY to ST
[st_primary_parameters, ~, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, original_parameteters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(primary_parameters_type_string, st_primary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results
assert(isequal(round(XY_parameters2{1},4),round(original_parameteters,4)));

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

% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(length(XY_parameters)==length(secondary_parameters));
for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(size(XY_parameters{ith_parameter}),size(secondary_parameters{ith_parameter})));
end

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters;
% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (original_parameteters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results

for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(round(XY_parameters2{ith_parameter},4),round(original_parameteters{ith_parameter},4)));
end


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

% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(length(XY_parameters)==length(secondary_parameters));
for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(size(XY_parameters{ith_parameter}),size(secondary_parameters{ith_parameter})));
end

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters;
% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (original_parameteters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results

for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(round(XY_parameters2{ith_parameter},4),round(original_parameteters{ith_parameter},4)));
end

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

% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check size of results
assert(iscell(XY_parameters));
assert(length(XY_parameters)==length(secondary_parameters));
for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(size(XY_parameters{ith_parameter}),size(secondary_parameters{ith_parameter})));
end

%%%
% Show that, now that parameters are fixed, can call the transform to ST
% and then back to XY and get the same results.

original_parameteters = XY_parameters;
% Call the function to convert from XY to ST
[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_primary_parameter_is_flipped] = ...
fcn_geometry_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (original_parameteters), (-1));

% Call the function to convert from ST back to XY
[XY_parameters2] = ...
fcn_geometry_orientGeometrySt2XY(secondary_parameters_type_strings, st_secondary_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, (fig_num));

% Check results

for ith_parameter = 1:length(secondary_parameters)
    assert(isequal(round(XY_parameters2{ith_parameter},4),round(original_parameteters{ith_parameter},4)));
end

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_orientGeometryXY2St(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end