%% script_test_fcn_geometry_alignArcArcC2Optimized
% Exercises the function: fcn_geometry_alignArcArcC2Optimized

% Revision history:
% 2024_07_26 S. Brennan
% -- wrote the code from fcn_geometry_alignArcArc

close all;


%% check input orientation corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Basic test 1.1 - checking plot inputs of arcs that are correctly oriented
fig_num = 11;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [0 2+0.2];
arc2_radius               = 2;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

% Put arcs into a fit sequence
fitSequence_bestFitType{1} = 'arc';
fitSequence_bestFitType{2} = 'arc';

fitSequence_parameters{1}  = arc1_parameters;
fitSequence_parameters{2}  = arc2_parameters;

% Join these using a spiral
continuity_level = 2;
[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

% Generate XY data
XY_data = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters,[], [], (fig_num));
plot(XY_data(:,1),XY_data(:,2),'k.');



% sgtitle('Checking that arc1 is joined to arc2: C1 continuous');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    3.0000   -3.1416   -1.5708         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.2000    3.0000    3.0000   -1.5708         0         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[1.0000         0         0         0         0    0.2000]));

%% Basic test 1.2 - checking plot inputs of arcs, arc1 is incorrectly oriented
fig_num = 12;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_end           = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [0.2 3];
arc2_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C1 continuous, arc 1 in bad orientation');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    3.0000   -3.1416   -1.5708         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.2000    3.0000    3.0000   -1.5708         0         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[1.0000         0         0         0         0    0.2000]));

%% Basic test 1.3 - checking plot inputs of arcs, arc2 is incorrectly oriented
fig_num = 13;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [0.2 3];
arc2_radius               = 3;
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_start         = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C1 continuous, arc 2 in bad orientation');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    3.0000   -3.1416   -1.5708         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.2000    3.0000    3.0000   -1.5708         0         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[1.0000         0         0         0         0    0.2000]));


%% Basic test 1.4 - checking plot inputs of arcs, arc1 and arc2 is incorrectly oriented
fig_num = 14;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_end           = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [0.2 3];
arc2_radius               = 3;
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_start         = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C1 continuous, arc 1 and arc2 in bad orientation');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    3.0000   -3.1416   -1.5708         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.2000    3.0000    3.0000   -1.5708         0         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[1.0000         0         0         0         0    0.2000]));


%% check conversions into St coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _______                                              _
%  / ____|__   __|                                            (_)
% | (___    | |     ______    ___ ___  _ ____   _____ _ __ ___ _  ___  _ __
%  \___ \   | |    |______|  / __/ _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
%  ____) |  | |             | (_| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
% |_____/   |_|              \___\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=ST%20-%20conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following section checks whether the ST conversion sub-codes are
% working

%% Basic test 2.1 - checking the + to + cross product combination
fig_num = 21;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.3;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 90*pi/180) sin(arc1_end_angle + 90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to + ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0    3.0000    1.0000   -2.0944   -0.7980         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.5830    3.0098    0.6000   -0.7980    1.0472         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.7160    0.6981    0.6981    2.2840         0    0.4243]));


%% Basic test 2.20 - checking the + to - cross product combination, no intersection
fig_num = 220;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -   90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -    0*pi/180) sin(arc1_end_angle -   0*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_radius               = 0.8;
offset_s = 0.3 ;
offset_t = -0.1 - 2*arc2_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle + 180*pi/180) sin(arc1_end_angle + 180*pi/180)];
arc2_vector_end           = [cos(arc1_end_angle +  90*pi/180) sin(arc1_end_angle +  90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to - ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    3.0000    1.0000   -2.0944   -0.7273         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 1.7954    2.3098    0.8000    2.4143    1.0472         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6649    0.7469    0.7469    2.3351         0    0.6782]));

%% Basic test 2.21 - checking the + to - cross product combination, intersection
fig_num = 221;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -   90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -    0*pi/180) sin(arc1_end_angle -   0*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_radius               = 0.8;
offset_s = 0.3 ;
offset_t = 0.1 - 2*arc2_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle + 180*pi/180) sin(arc1_end_angle + 180*pi/180)];
arc2_vector_end           = [cos(arc1_end_angle +  90*pi/180) sin(arc1_end_angle +  90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to - ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    3.0000    1.0000   -2.0944   -0.3489         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 1.6915    2.3846    0.8000    2.7927    1.0472         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.3419    0.9397    0.9397    2.6581         0         0]));

%% Basic test 2.30 - checking the - to + cross product combination, no intersection
fig_num = 230;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -  90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -  180*pi/180) sin(arc1_end_angle - 180*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles               = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2)    = arc1_center_xy;
arc1_parameters(1,3)      = arc1_radius;
arc1_parameters(1,4:5)    = arc1_angles;
arc1_parameters(1,6)      = arc1_is_circle;
arc1_parameters(1,7)      = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.3;
offset_t = 0.1 + 2*arc1_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.8;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 90*pi/180) sin(arc1_end_angle + 90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0    3.0000    1.0000   -2.0944    2.8217         0         0]));
assert(isequal(round(revised_arc2_parameters,4),[ -1.4954    4.2098    0.8000   -0.3199    1.0472         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.3144    0.9493   -0.9493    3.3144         0    0.6782]));


%% Basic test 2.311 - checking the - to + cross product combination, intersecting circles, no intersecting arcs
fig_num = 2311;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -  90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -  180*pi/180) sin(arc1_end_angle - 180*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles               = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2)    = arc1_center_xy;
arc1_parameters(1,3)      = arc1_radius;
arc1_parameters(1,4:5)    = arc1_angles;
arc1_parameters(1,6)      = arc1_is_circle;
arc1_parameters(1,7)      = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.3;
offset_t = -0.1 + 2*arc1_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.8;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 90*pi/180) sin(arc1_end_angle + 90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    3.0000    1.0000   -2.0944    2.4433         0         0]));
assert(isequal(round(revised_arc2_parameters,4),[ -1.3787    4.1572    0.8000   -0.6983    1.0472         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6429    0.7660   -0.7660    3.6429         0         0]));


%% Basic test 2.312 - checking the - to + cross product combination, intersecting circles, no intersecting arcs
fig_num = 2312;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -  90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -  180*pi/180) sin(arc1_end_angle - 180*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles               = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2)    = arc1_center_xy;
arc1_parameters(1,3)      = arc1_radius;
arc1_parameters(1,4:5)    = arc1_angles;
arc1_parameters(1,6)      = arc1_is_circle;
arc1_parameters(1,7)      = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.3;
offset_t = -0.1 + 2*arc1_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.8;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle - 40*pi/180) sin(arc1_end_angle - 40*pi/180)];
arc2_vector_end           = [cos(arc1_end_angle + 90*pi/180) sin(arc1_end_angle + 90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    3.0000    1.0000   -2.0944    2.4433         0         0]));
assert(isequal(round(revised_arc2_parameters,4),[ -1.3787    4.1572    0.8000   -0.6983    1.0472         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6429    0.7660   -0.7660    3.6429         0         0]));


%% Basic test 2.313 - checking the - to + cross product combination, intersecting circles, two intersecting arcs
fig_num = 233;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -  90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -  180*pi/180  - 40*pi/180) sin(arc1_end_angle - 180*pi/180  - 40*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles               = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2)    = arc1_center_xy;
arc1_parameters(1,3)      = arc1_radius;
arc1_parameters(1,4:5)    = arc1_angles;
arc1_parameters(1,6)      = arc1_is_circle;
arc1_parameters(1,7)      = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.3;
offset_t = -0.1 + 2*arc1_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.8;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle - 40*pi/180) sin(arc1_end_angle - 40*pi/180)];
arc2_vector_end           = [cos(arc1_end_angle + 90*pi/180) sin(arc1_end_angle + 90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    3.0000    1.0000   -2.0944    2.4433         0         0]));
assert(isequal(round(revised_arc2_parameters,4),[ -1.3787    4.1572    0.8000   -0.6983    1.0472         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6429    0.7660   -0.7660    3.6429         0         0]));



%% Basic test 2.4 - checking the - to - cross product combination
fig_num = 24;
figure(fig_num); clf;

tolerance = 0.4; % meters

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -   90*pi/180) sin(arc1_end_angle -  90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle -  180*pi/180) sin(arc1_end_angle - 180*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles               = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2)    = arc1_center_xy;
arc1_parameters(1,3)      = arc1_radius;
arc1_parameters(1,4:5)    = arc1_angles;
arc1_parameters(1,6)      = arc1_is_circle;
arc1_parameters(1,7)      = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_radius               = 0.8;
offset_s = 0.3 ;
offset_t = 0.1; %- 2*arc2_radius;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle + 180*pi/180) sin(arc1_end_angle + 180*pi/180)];
arc2_vector_end           = [cos(arc1_end_angle +  90*pi/180) sin(arc1_end_angle +  90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to - ');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944    1.6115         0         0]));
assert(isequal(round(revised_arc2_parameters,4),[0.2366    3.2098    0.8000    1.6115    1.0472         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.9992    0.0407   -0.0407    3.9992         0    0.2449]));



%% check C0 intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _               _    _                _____ ___
%  / ____| |             | |  (_)              / ____/ _ \
% | |    | |__   ___  ___| | ___ _ __   __ _  | |   | | | |
% | |    | '_ \ / _ \/ __| |/ / | '_ \ / _` | | |   | | | |
% | |____| | | |  __/ (__|   <| | | | | (_| | | |___| |_| |
%  \_____|_| |_|\___|\___|_|\_\_|_| |_|\__, |  \_____\___/
%                                       __/ |
%  _____       _                       |___/  _
% |_   _|     | |                        | | (_)
%   | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Checking%20C0%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3


%% Basic test 3.11 - checking the + to + cross product combination, large to small, feasible, no intersection
fig_num = 311;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = 0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;

% revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters

[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3464    2.8000    0.6000   -0.5236    1.5708         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));



%% Basic test 3.12 - checking the + to + cross product combination, large to small with shift blocked, no intersection
% Setting tolerance to empty makes shift blocked
fig_num = 312;
figure(fig_num); clf;

tolerance = []; % meters

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 3.131 - checking the + to + cross product combination, large to small with shift allowed, no intersection
fig_num = 3131;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3464    2.8000    0.6000   -0.5236    1.5708         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 3.132 - checking the + to + cross product combination, large to small with shift allowed, one intersection
fig_num = 3131;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = 30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_start_angle          = -30*pi/180;
arc2_vector_start         = [cos(arc2_start_angle) sin(arc2_start_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -1.0472    0.0501         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.4330    3.2500    0.6000   -0.3396    2.6180         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 3.133 - checking the + to + cross product combination, large to small with shift allowed, two intersections
fig_num = 3133;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = 30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle +50*pi/180) sin(arc1_end_angle +50*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_start_angle          = -30*pi/180;
arc2_vector_start         = [cos(arc2_start_angle) sin(arc2_start_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -1.0472    0.0501         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.4330    3.2500    0.6000   -0.3396    2.6180         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));



%% Basic test 3.21 - checking the + to + cross product combination, small to large, feasible
fig_num = 321;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.7321    1.6000    3.0000   -0.5236    1.5708         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.22 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 322;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.23 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 323;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.7321    1.6000    3.0000   -0.5236    1.5708         0    1.0000]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));



%% Basic test 3.31 - checking the + to - cross product combination, large to small, feasible
fig_num = 331;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    3.0000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 1.3856    2.2000    0.6000    2.6180   -1.5708         0         0]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.32 - checking the + to - cross product combination, large to small with shift blocked
% Setting tolerance to empty makes shift blocked
fig_num = 332;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.33 - checking the + to - cross product combination, large to small with shift allowed
fig_num = 333;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[1.3856    2.2000    0.6000    2.6180   -1.5708         0         0]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));




%% Basic test 3.41 - checking the + to + cross product combination, small to large, feasible
fig_num = 341;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 3.4641   -1.4000    3.0000    2.6180   -1.5708         0         0]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.42 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 342;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.43 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 343;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 0;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 3.4641   -1.4000    3.0000    2.6180   -1.5708         0         0]));
assert(isempty(revised_intermediate_geometry_join_type));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% check C1 intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _               _    _                _____ __
%  / ____| |             | |  (_)              / ____/_ |
% | |    | |__   ___  ___| | ___ _ __   __ _  | |     | |
% | |    | '_ \ / _ \/ __| |/ / | '_ \ / _` | | |     | |
% | |____| | | |  __/ (__|   <| | | | | (_| | | |____ | |
%  \_____|_| |_|\___|\___|_|\_\_|_| |_|\__, |  \_____||_|
%                                       __/ |
%  _____       _                       |___/  _
% |_   _|     | |                        | | (_)
%   | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Checking%20C1%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4

%% Basic test 4.11 - checking the + to + cross product combination, large to small, feasible
fig_num = 411;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = 0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;

% revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters

[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.3817         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3712    2.8510    0.6000   -0.3817    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.3725    0.9280    0.9280    2.6275         0         0]));

%% Basic test 4.12 - checking the + to + cross product combination, large to small with shift blocked
% Setting tolerance to empty makes shift blocked
fig_num = 412;
figure(fig_num); clf;

tolerance = []; % meters

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 4.13 - checking the + to + cross product combination, large to small with shift allowed
fig_num = 413;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -1.1671         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.4330    2.7500    0.6000   -1.1671    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.9196    0.3928    0.3928    2.0804         0    0.3000]));




%% Basic test 4.21 - checking the + to + cross product combination, small to large, feasible
fig_num = 421;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    0.6000    1.0000   -2.0944   -0.5499         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.7051    1.6452    3.0000   -0.5499    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.5226    0.8526    0.8526    0.0774         0         0]));

%% Basic test 4.22 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 422;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.23 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 423;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.7321    1.6000    3.0000   -0.5236    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.5000    0.8660    0.8660    0.1000         0    0.0000]));



%% Basic test 4.31 - checking the + to - cross product combination, large to small, feasible
fig_num = 431;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[-0.0000    3.0000    1.0000   -2.0944   -0.7419         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[1.4539    2.2183    0.6000    2.3997   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6757    0.7372    0.7372    2.3243         0    0.4062]));

%% Basic test 4.32 - checking the + to - cross product combination, large to small with shift blocked
% Setting tolerance to empty makes shift blocked
fig_num = 432;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.33 - checking the + to - cross product combination, large to small with shift allowed
fig_num = 433;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[1.3856    2.2000    0.6000    2.6180   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.5000    0.8660    0.8660    2.5000         0         0]));




%% Basic test 4.41 - checking the + to + cross product combination, small to large, feasible
fig_num = 441;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0    0.6000    1.0000   -2.0944   -0.7330         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[3.5757   -1.4067    3.0000    2.4085   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6691    0.7431    0.7431   -0.0691         0    0.9014]));

%% Basic test 4.42 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 442;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.43 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 443;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 1;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5236         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 3.4641   -1.4000    3.0000    2.6180   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[ 0.5000    0.8660    0.8660    0.1000         0         0]));




%% check C2 intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _               _    _                _____ ___
%  / ____| |             | |  (_)              / ____|__ \
% | |    | |__   ___  ___| | ___ _ __   __ _  | |       ) |
% | |    | '_ \ / _ \/ __| |/ / | '_ \ / _` | | |      / /
% | |____| | | |  __/ (__|   <| | | | | (_| | | |____ / /_
%  \_____|_| |_|\___|\___|_|\_\_|_| |_|\__, |  \_____|____|
%                                       __/ |
%  _____       _                       |___/  _
% |_   _|     | |                        | | (_)
%   | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Checking%20C2%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5

%% Basic test 5.11 - checking the + to + cross product combination, large to small, feasible
fig_num = 511;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = 0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    1.0000   -2.0944   -1.0311         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3281    2.8683    0.6000    0.7311    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.5139    2.1421    0.5397    1.3217    1.0000    1.6667]));

%% Basic test 5.12 - checking the + to + cross product combination, large to small with shift blocked
% Setting tolerance to empty makes shift blocked
fig_num = 512;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.13 - checking the + to + cross product combination, large to small with shift allowed
fig_num = 513;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.6185         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3455    2.8005    0.6000   -0.3654    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[ 0.8148    2.4202    0.9523    0.1898    1.0000    1.6667]));




%% Basic test 5.21 - checking the + to + cross product combination, small to large, feasible
fig_num = 521;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -1.5165         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.6204    1.5933    3.0000   -0.2366    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.0542   -0.3985    0.0543    1.9199    1.0000    0.3333]));

%% Basic test 5.22 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 522;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.23 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 523;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start - arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180) sin(arc1_end_angle + 120*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.6185         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[-1.7312    1.5995    3.0000   -0.4920    1.5708         0    1.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.8148    0.0202    0.9523    0.1897    1.0000    0.3333]));



%% Basic test 5.31 - checking the + to - cross product combination, large to small, feasible
fig_num = 531;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.05;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.8289         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[1.4539    2.2183    0.6000    2.0852   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.6757    2.2628    0.7419    0.6825    1.0000   -1.6667]));

%% Basic test 5.32 - checking the + to - cross product combination, large to small with shift blocked
% Setting tolerance to empty makes shift blocked
fig_num = 532;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.33 - checking the + to - cross product combination, large to small with shift allowed
fig_num = 533;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 0.6;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.5710         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[1.3865    2.1995    0.6000    2.5389   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.8413    2.4595    0.9998    0.0949    1.0000   -1.6667]));




%% Basic test 5.41 - checking the + to + cross product combination, small to large, feasible
fig_num = 541;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0.05;
offset_t = -0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -1.1859         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[3.5757   -1.4067    3.0000    2.4083   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.3754   -0.3269    0.3849    1.3580    1.0000   -0.3333]));

%% Basic test 5.42 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 542;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.43 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 543;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 0.6];
arc1_radius               = 1;
arc1_end_angle            = -30*pi/180;
arc1_vector_start         = [cos(arc1_end_angle -90*pi/180) sin(arc1_end_angle -90*pi/180)];
arc1_vector_end           = [cos(arc1_end_angle) sin(arc1_end_angle)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
offset_s = 0;
offset_t = 0.1;

start_unit_radial_vector  = [cos(arc1_end_angle) sin(arc1_end_angle)];
start_unit_tangent_vector = start_unit_radial_vector*[0 1; -1 0];

end_point_arc1 = arc1_center_xy + start_unit_radial_vector*arc1_radius;
offset_point_arc2_start = end_point_arc1 + start_unit_tangent_vector*offset_s -start_unit_radial_vector*offset_t;
arc2_radius               = 3;
arc2_center_xy            = offset_point_arc2_start + arc2_radius*start_unit_radial_vector;

arc2_vector_start         = [cos(arc1_end_angle+pi) sin(arc1_end_angle+pi)];
arc2_vector_end           = [cos(arc1_end_angle + 120*pi/180+pi) sin(arc1_end_angle + 120*pi/180+pi)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcArcC2Optimized(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[ 0.0000    0.6000    1.0000   -2.0944   -0.5907         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[ 3.4650   -1.4005    3.0000    2.5956   -1.5708         0         0]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.8306    0.0431    0.9801    0.1342    1.0000   -0.3333]));






%% Fail conditions
if 1==0
    % FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_alignArcArcC2Optimized(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end