%% script_test_fcn_geometry_flipGeom
% Exercises the function: fcn_geometry_flipGeom

% Revision history:
% 2024_05_15 - Sean Brennan
% -- wrote the code
% 2024_06_19 - Sean Brennan
% -- changed parameter format to new style:
%            'spiral' - 
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
% 2024_06_19 - Sean Brennan
% -- changed parameter format for line to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
% 2024_06_19 - Sean Brennan
% -- changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]

% TO-DO
% (none)

close all

%% Basic Test: circle - no change
fig_num = 1;
figure(fig_num); clf;

% Fill in circle 1
circle1_center_xy            = [-3 3];
circle1_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

geomType = 'circle';
geomParameters = circle1_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(isequal(round(geomParameters_flipped,4), geomParameters));

%% Basic Test: arc - flips start and end angles, and clockwise or CCW flag
fig_num = 2;
figure(fig_num); clf;

% Fill in arc 2
arc2_center_xy            = [3 0];
arc2_radius               = 2;
arc2_vector_start         = [cos(-135*pi/180) sin(-135*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

geomType = 'arc';
geomParameters = arc2_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(isequal(round(geomParameters_flipped,4), round([arc2_center_xy arc2_radius fliplr(arc2_angles') arc2_is_circle ~arc2_is_counter_clockwise],4)));

%% Basic Test: line - flips direction
fig_num = 3;
figure(fig_num); clf;

% Fill in line
line_unit_tangent_vector = [0 1];
line_base_point_xy       = [-4 3];

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3)   = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

geomType = 'line';
geomParameters = line_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(isequal(round(geomParameters_flipped,4), round([line_base_point_xy line_parameters(1,3)+pi],4)));

%% Basic Test: segment - flips start and end stations, flips direction, but otherwise the same
fig_num = 4;
figure(fig_num); clf;

% Fill in line
segment_unit_tangent_vector = [0 1];
segment_base_point_xy       = [-4 3];
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3)   = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

geomType = 'segment';
geomParameters = segment_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(isequal(round(geomParameters_flipped,4), round([segment_base_point_xy+segment_unit_tangent_vector*segment_length segment_parameters(1,3)+pi segment_parameters(1,4)],4)));

%% Basic Test: spiral - flips start and end stations, flips direction, but otherwise the same
fig_num = 5;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear spiral_parameters
spiral_parameters = [2.2403   -1.1202   -0.9002    0.5645    1.0000         0];

geomType = 'spiral';
geomParameters = spiral_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(isequal(round(geomParameters_flipped,4), round([2.6662   -1.4877    2.5236    0.5645         0   -1.0000],4)));

%% Basic Test: 'none' - does nothing
fig_num = 6;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear none_parameters
none_parameters = [nan nan];

geomType = 'none';
geomParameters = none_parameters;

geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_flipped)));
assert(all(isnan(none_parameters)));




