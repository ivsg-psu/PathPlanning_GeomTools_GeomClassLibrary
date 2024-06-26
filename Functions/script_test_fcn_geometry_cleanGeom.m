%% script_test_fcn_geometry_cleanGeom
% Exercises the function: fcn_geometry_cleanGeom

% Revision history:
% 2024_06_25 - Sean Brennan
% -- wrote the code


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
geomParameters_good = circle1_parameters;

bad_circle1_parameters = circle1_parameters;

geomType = 'circle';
geomParameters = bad_circle1_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(isequal(round(geomParameters_good,4),  round(geomParameters_cleaned,4)));

%% Basic Test: arc - bad start and end angles
fig_num = 2;
figure(fig_num); clf;

% Fill in arc 2
arc2_center_xy            = [3 0];
arc2_radius               = 2;
arc2_vector_start         = [cos(-135*pi/180) sin(-135*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = mod([atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));],2*pi);


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;
geomParameters_good = arc2_parameters;

% Make bad arc
bad_arc2_parameters = arc2_parameters;
bad_arc2_parameters(1,4:5) = arc2_parameters(1,4:5)-2*pi;

geomType = 'arc';
geomParameters = bad_arc2_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(isequal(round(geomParameters_good,4),  round(geomParameters_cleaned,4)));

%% Basic Test: line - bad angle
fig_num = 3;
figure(fig_num); clf;

% Fill in line
line_unit_tangent_vector = [0 1];
line_base_point_xy       = [-4 3];

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3)   = mod(atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1)),2*pi);
geomParameters_good = line_parameters;

bad_line_parameters = line_parameters;
bad_line_parameters(1,3) = line_parameters(1,3)-2*pi;

geomType = 'line';
geomParameters = bad_line_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(isequal(round(geomParameters_good,4),  round(geomParameters_cleaned,4)));

%% Basic Test: segment - bad angle and length
fig_num = 4;
figure(fig_num); clf;

% Fill in line
segment_unit_tangent_vector = [0 1];
segment_base_point_xy       = [-4 3];
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3)   = mod(atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1)),2*pi);
segment_parameters(1,4)   = segment_length;
geomParameters_good = segment_parameters;

bad_segment_parameters(1,1:2) = segment_base_point_xy + segment_length*segment_unit_tangent_vector;
bad_segment_parameters(1,3)   = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1)) - 2*pi;
bad_segment_parameters(1,4)   = -1*segment_length;

geomType = 'segment';
geomParameters = bad_segment_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(isequal(round(geomParameters_good,4),  round(geomParameters_cleaned,4)));

%% Basic Test: spiral - flips start and end stations, flips direction, but otherwise the same
fig_num = 5;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear spiral_parameters
spiral_parameters = [2.2403   -1.1202   mod(-0.9002,2*pi)    0.5645    1.0000         0];
geomParameters_good = spiral_parameters;


bad_spiral_parameters = spiral_parameters;
spiral_parameters(1,3) = spiral_parameters(1,3) - 2*pi;

geomType = 'spiral';
geomParameters = bad_spiral_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(isequal(round(geomParameters_good,4),  round(geomParameters_cleaned,4)));

%% Basic Test: 'none' - does nothing
fig_num = 6;
figure(fig_num); clf;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear none_parameters
none_parameters = [nan nan];

geomType = 'none';
geomParameters = none_parameters;

geomParameters_cleaned = fcn_geometry_cleanGeom(geomType,  geomParameters, fig_num);

assert(isequal(size(geomParameters),size(geomParameters_cleaned)));
assert(all(isnan(none_parameters)));




