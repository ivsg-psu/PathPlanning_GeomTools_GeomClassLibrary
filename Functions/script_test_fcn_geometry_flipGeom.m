%% script_test_fcn_geometry_flipGeom
% Exercises the function: fcn_geometry_flipGeom
% Revision history:
% 2024_05_15 - Sean Brennan
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

%% Basic Test: line - flips start and end points, but otherwise the same
fig_num = 3;
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




