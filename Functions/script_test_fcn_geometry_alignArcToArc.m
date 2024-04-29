%% script_test_fcn_geometry_alignArcToArc
% Exercises the function: fcn_geometry_alignArcToArc
% Revision history:
% 2024_04_21
% -- wrote the code

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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C2 continuous');


% Check size of results
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_arc1_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));

% Check the arc results
fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));

% Check the spiral results
assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));


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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C2 continuous');

% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% % true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% % true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));
% 
% % Check the spiral results
% assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));



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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C2 continuous');

% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% % true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% % true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));
% 
% % Check the spiral results
% assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));


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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking that arc1 is joined to arc2: C2 continuous');

% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% % true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% % true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));
% 
% % Check the spiral results
% assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));


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
offset_t = 0.1;

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to + ');


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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to - ');

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to - ');

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');

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

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to + ');



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


continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are - to - ');

%%%%%%%%%%%%%%%%%%%%
% 
% %% Basic test 1.1 - an arc nearby the line segment joined with C0 continuity
% fig_num = 11;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% shift_error = [0 0.2];
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% 
% line_unit_tangent_vector = true_line_unit_tangent_vector;
% line_base_point_xy       = true_start_point_xy;
% line_s_start             = 0;
% line_s_end               = 1;
% 
% 
% true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;
% 
% arc_center_xy            = true_arc_center_xy+shift_error;
% arc_radius               = 1;
% arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
% arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% flag_arc_is_first = 0;
% continuity_level = 0;
% [revised_arc1_parameters, revised_arc2_parameters] = fcn_geometry_alignArcToArc(...
%     line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
% 
% title('Checking that arc is joined to the line: C0 continuous');
% 
% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal(fitted_line_params, true_line_params));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% % true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal([-0.3420    0.9397    1.0000    5.0615         0         0    1.0000], fitted_arc_params));
% 
% %% Basic test 1.2 - an arc intersecting the line segment joined with C0 continuity
% fig_num = 12;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% shift_error = [-0.8 -0.2];
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% 
% line_unit_tangent_vector = true_line_unit_tangent_vector;
% line_base_point_xy       = true_start_point_xy;
% line_s_start             = 0;
% line_s_end               = 1;
% 
% 
% true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;
% 
% arc_center_xy            = true_arc_center_xy+shift_error;
% arc_radius               = 1;
% arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
% arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% flag_arc_is_first = 0;
% continuity_level = 0;
% [revised_arc1_parameters, revised_arc2_parameters] = fcn_geometry_alignArcToArc(...
%     line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
% 
% title('Checking that arc is joined to the line: C0 continuous');
% 
% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% % true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal([1     0    -1     0     0     1 ], fitted_line_params));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% % true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal([-0.3420    0.9397    1.0000    5.0615         0         0    1.0000], fitted_arc_params));
% 
% %% ADVANCED test 1.9 - systematic joining line to arc, C0 continuity
% fig_num = 19;
% 
% % Bit 1 (top bit)            : line precedes arc, 
% % Bit 2 (2nd from top bit)   : line oriented, 
% % Bit 3                      : arc oriented, 
% % Bit 4                      : ends aligned, 
% % Bit 5 (lowest bit)         : arc is positive
% NtotalTests = (2^5);
% for ith_test = 1:NtotalTests
% 
%     numDigits = 5;
%     binary_string = dec2bin(ith_test-1,numDigits);
% 
%     figure(fig_num); clf;
%     title_string = sprintf('Test %.0d of %.0d, %s: ',ith_test, NtotalTests, binary_string);
% 
% 
%     % Top bit - is the line first? 0 is yes, 1 is no
%     if strcmp(binary_string(1),'0')
%         title_string = cat(2,title_string,'line precedes arc,');
% 
%         true_line_unit_tangent_vector = [1 0];
%         true_start_point_xy = [-1 0];
% 
%         line_unit_tangent_vector = [1 0];
%         line_base_point_xy       = [-1 0];
%         line_s_start             = 0;
%         line_s_end               = 1;
% 
%         arc_center_xy            = [0 1];
%         arc_radius               = 1;
%         arc_vector_start         = [ 0 -1];
%         arc_vector_end           = [ 1  0];
%         arc_is_circle            = 0;
%         arc_is_counter_clockwise = 1;
% 
%         true_arc_center_xy  = [0 1];
%         true_arc_is_counter_clockwise = 1;
%         true_arc_angles     = [270 360]*pi/180;
% 
%         flag_arc_is_first = 0;
%     else
%         title_string = cat(2,title_string,'arc precedes line,');
% 
%         true_line_unit_tangent_vector = [1 0];
%         true_start_point_xy           = [0 0];
% 
%         line_unit_tangent_vector = [1 0];
%         line_base_point_xy       = [0 0];
%         line_s_start             = 0;
%         line_s_end               = 1;
% 
%         arc_center_xy            = [0 1];
%         arc_radius               = 1;
%         arc_vector_start         = [-1  0];
%         arc_vector_end           = [ 0 -1];  
%         arc_is_circle            = 0;
%         arc_is_counter_clockwise = 1;
% 
%         true_arc_center_xy  = [0 1];
%         true_arc_is_counter_clockwise = 1;
%         true_arc_angles          = [180 270]*pi/180;
%         flag_arc_is_first = 1;
%     end
% 
%     % Next from top bit - is the line oriented correctly? 0 is yes, 1 is no
%     if strcmp(binary_string(2),'0')
%         title_string = cat(2,title_string,'line oriented,');
%     else
%         title_string = cat(2,title_string,'line misoriented,');
%         % Move the start point of the line to the end
%         line_base_point_xy = line_base_point_xy + line_s_end*line_unit_tangent_vector;
% 
%         % Change the line's orientation
%         line_unit_tangent_vector = -line_unit_tangent_vector;
%     end
% 
%     % Next from top bit - is the arc oriented correctly? 0 is yes, 1 is no
%     if strcmp(binary_string(3),'0')
%         title_string = cat(2,title_string,'arc oriented,');
%     else
%         title_string = cat(2,title_string,'arc misoriented,');
%         temp = arc_vector_start;
%         arc_vector_start = arc_vector_end;
%         arc_vector_end   = temp;
%         arc_is_counter_clockwise = ~arc_is_counter_clockwise;
%     end
% 
%     % Next from top bit - are the ends aligned? 0 is yes, 1 is no
%     if strcmp(binary_string(4),'0')
%         title_string = cat(2,title_string,'ends aligned,');
%     else
%         title_string = cat(2,title_string,'ends misaligned,');
%         if strcmp(binary_string(1),'0')
%             % Line is first, misalign the arc
%             arc_center_xy = arc_center_xy+[0 0.2];
%         else
%             % Arc is first, misalign the line
%             line_base_point_xy = line_base_point_xy + [0 0.2];
%         end
%     end
% 
%     % Bottom bit - is the arc oriented to the left? 0 is yes, 1 is no
%     if strcmp(binary_string(5),'0')
%         title_string = cat(2,title_string,'arc left,');
%     else
%         title_string = cat(2,title_string,'arc right,');
% 
%         % Flip the y values
%         arc_vector_start(2) = -arc_vector_start(2);
%         arc_vector_end(2)   = -arc_vector_end(2);
%         arc_is_counter_clockwise = ~arc_is_counter_clockwise;
%         arc_center_xy = -1*arc_center_xy;
% 
%         true_arc_center_xy = -1*true_arc_center_xy;        
%         true_arc_is_counter_clockwise = 0;
%         true_arc_angles = 2*pi-true_arc_angles;
%     end
% 
%     arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
%     % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
%     line_parameters(1,1:2) = line_unit_tangent_vector;
%     line_parameters(1,3:4) = line_base_point_xy;
%     line_parameters(1,5)   = line_s_start;
%     line_parameters(1,6)   = line_s_end;
% 
%     % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
%     arc_parameters(1,1:2) = arc_center_xy;
%     arc_parameters(1,3)   = arc_radius;
%     arc_parameters(1,4:5) = arc_angles;
%     arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
%     tolerance = 0.4; % Use default
% 
% 
%     continuity_level = 0;
%     [revised_arc1_parameters, revised_arc2_parameters] = ...
%         fcn_geometry_alignArcToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
%     sgtitle(title_string);
%     pause(0.01);
% 
%     % Check size of results
%     assert(isequal(size(revised_arc1_parameters),[1 7]));
%     assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
%     % Check that line results
%     fitted_line_params = round(revised_arc1_parameters,4);
%     true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
%     assert(isequal(fitted_line_params, true_line_params));
% 
%     % Check the arc results
%     fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
%     true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
%     assert(isequal(true_arc_params, fitted_arc_params));
% end
% 
% 
% 
% %% Basic test 2.1 - an arc nearby the line segment joined with C1 continuity
% fig_num = 21;
% figure(fig_num); clf;
% 
% tolerance = 0.4; % meters
% shift_error = [0 0.2];
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% 
% line_unit_tangent_vector = true_line_unit_tangent_vector;
% line_base_point_xy       = true_start_point_xy;
% line_s_start             = 0;
% line_s_end               = 1;
% 
% 
% true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;
% 
% arc_center_xy            = true_arc_center_xy+shift_error;
% arc_radius               = 1;
% arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
% arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% flag_arc_is_first = 0;
% continuity_level = 1;
% [revised_arc1_parameters, revised_arc2_parameters] = fcn_geometry_alignArcToArc(...
%     line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
% 
% title('Checking that arc is joined to the line: C1 continuous');
% 
% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal(fitted_line_params, true_line_params));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal(true_arc_params, fitted_arc_params));
% 
% %% Basic test 2.2 - an arc intersecting the line segment joined with C1 continuity
% fig_num = 22;
% figure(fig_num); clf;
% 
% tolerance = 1; % meters
% shift_error = [-0.8 -0.2];
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% 
% line_unit_tangent_vector = true_line_unit_tangent_vector;
% line_base_point_xy       = true_start_point_xy;
% line_s_start             = 0;
% line_s_end               = 1;
% 
% 
% true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;
% 
% arc_center_xy            = true_arc_center_xy+shift_error;
% arc_radius               = 1;
% arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
% arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% flag_arc_is_first = 0;
% continuity_level = 1;
% [revised_arc1_parameters, revised_arc2_parameters] = fcn_geometry_alignArcToArc(...
%     line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
% 
% title('Checking that arc is joined to the line: C0 continuous');
% 
% 
% % Check size of results
% assert(isequal(size(revised_arc1_parameters),[1 7]));
% assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
% % Check that line results
% fitted_line_params = round(revised_arc1_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
% assert(isequal(round(true_line_params,4), fitted_line_params));
% 
% % Check the arc results
% fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
% assert(isequal(round(true_arc_params,4), fitted_arc_params));
% 
% 
% %% ADVANCED test 2.9 - joining line to arc, C1 continuity
% fig_num = 29;
% 
% % Bit 1 (top bit)            : line precedes arc, 
% % Bit 2 (2nd from top bit)   : line oriented, 
% % Bit 3                      : arc oriented, 
% % Bit 4                      : ends aligned, 
% % Bit 5 (lowest bit)         : arc is positive
% NtotalTests = (2^5);
% for ith_test = 1:NtotalTests
% 
%     numDigits = 5;
%     binary_string = dec2bin(ith_test-1,numDigits);
% 
%     figure(fig_num); clf;
%     title_string = sprintf('Test %.0d of %.0d, %s: ',ith_test, NtotalTests, binary_string);
%     title(title_string)
% 
%     % Top bit - is the line first? 0 is yes, 1 is no
%     if strcmp(binary_string(1),'0')
%         title_string = cat(2,title_string,'line precedes arc,');
% 
%         true_line_unit_tangent_vector = [1 0];
%         true_start_point_xy = [-1 0];
% 
%         line_unit_tangent_vector = [1 0];
%         line_base_point_xy       = [-1 0];
%         line_s_start             = 0;
%         line_s_end               = 1;
% 
%         arc_center_xy            = [0 1];
%         arc_radius               = 1;
%         arc_vector_start         = [ 0 -1];
%         arc_vector_end           = [ 1  0];
%         arc_is_circle            = 0;
%         arc_is_counter_clockwise = 1;
% 
%         true_arc_center_xy  = [0 1];
%         true_arc_is_counter_clockwise = 1;
%         true_arc_angles     = [270 360]*pi/180;
% 
%         flag_arc_is_first = 0;
%     else
%         title_string = cat(2,title_string,'arc precedes line,');
% 
%         true_line_unit_tangent_vector = [1 0];
%         true_start_point_xy           = [0 0];
% 
%         line_unit_tangent_vector = [1 0];
%         line_base_point_xy       = [0 0];
%         line_s_start             = 0;
%         line_s_end               = 1;
% 
%         arc_center_xy            = [0 1];
%         arc_radius               = 1;
%         arc_vector_start         = [-1  0];
%         arc_vector_end           = [ 0 -1];  
%         arc_is_circle            = 0;
%         arc_is_counter_clockwise = 1;
% 
%         true_arc_center_xy  = [0 1];
%         true_arc_is_counter_clockwise = 1;
%         true_arc_angles          = [180 270]*pi/180;
%         flag_arc_is_first = 1;
%     end
% 
%     % Next from top bit - is the line oriented correctly? 0 is yes, 1 is no
%     if strcmp(binary_string(2),'0')
%         title_string = cat(2,title_string,'line oriented,');
%     else
%         title_string = cat(2,title_string,'line misoriented,');
%         % Move the start point of the line to the end
%         line_base_point_xy = line_base_point_xy + line_s_end*line_unit_tangent_vector;
% 
%         % Change the line's orientation
%         line_unit_tangent_vector = -line_unit_tangent_vector;
%     end
% 
%     % Next from top bit - is the arc oriented correctly? 0 is yes, 1 is no
%     if strcmp(binary_string(3),'0')
%         title_string = cat(2,title_string,'arc oriented,');
%     else
%         title_string = cat(2,title_string,'arc misoriented,');
%         temp = arc_vector_start;
%         arc_vector_start = arc_vector_end;
%         arc_vector_end   = temp;
%         arc_is_counter_clockwise = ~arc_is_counter_clockwise;
%     end
% 
%     % Next from top bit - are the ends aligned? 0 is yes, 1 is no
%     if strcmp(binary_string(4),'0')
%         title_string = cat(2,title_string,'ends aligned,');
%     else
%         title_string = cat(2,title_string,'ends misaligned,');
%         if strcmp(binary_string(1),'0')
%             % Line is first, misalign the arc
%             arc_center_xy = arc_center_xy+[0 0.2];
%         else
%             % Arc is first, misalign the line
%             line_base_point_xy = line_base_point_xy + [0 0.2];
%         end
%     end
% 
%     % Bottom bit - is the arc oriented to the left? 0 is yes, 1 is no
%     if strcmp(binary_string(5),'0')
%         title_string = cat(2,title_string,'arc left,');
%     else
%         title_string = cat(2,title_string,'arc right,');
% 
%         % Flip the y values
%         arc_vector_start(2) = -arc_vector_start(2);
%         arc_vector_end(2)   = -arc_vector_end(2);
%         arc_is_counter_clockwise = ~arc_is_counter_clockwise;
%         arc_center_xy = -1*arc_center_xy;
% 
%         true_arc_center_xy = -1*true_arc_center_xy;        
%         true_arc_is_counter_clockwise = 0;
%         true_arc_angles = 2*pi-true_arc_angles;
%     end
% 
%     arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
%     % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
%     line_parameters(1,1:2) = line_unit_tangent_vector;
%     line_parameters(1,3:4) = line_base_point_xy;
%     line_parameters(1,5)   = line_s_start;
%     line_parameters(1,6)   = line_s_end;
% 
%     % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
%     arc_parameters(1,1:2) = arc_center_xy;
%     arc_parameters(1,3)   = arc_radius;
%     arc_parameters(1,4:5) = arc_angles;
%     arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
%     tolerance = 0.4; % Use default
% 
% 
%     continuity_level = 1;
%     [revised_arc1_parameters, revised_arc2_parameters] = ...
%         fcn_geometry_alignArcToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
%     sgtitle(title_string);
%     pause(0.01);
% 
%     % Check size of results
%     assert(isequal(size(revised_arc1_parameters),[1 7]));
%     assert(isequal(size(revised_arc2_parameters),[1 7]));
% 
%     % Check that line results
%     fitted_line_params = round(revised_arc1_parameters,4);
%     true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
%     assert(isequal(fitted_line_params, true_line_params));
% 
%     % Check the arc results
%     fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
%     true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
%     assert(isequal(true_arc_params, fitted_arc_params));
% end
% 
% %% Test that threshold will not cause fitting if error is too large
% fig_num = 2;
% figure(fig_num); clf;
% 
% tolerance = 0.01; % Use very low tolerance to force fit to NOT occur
% shift_error = [0 0.2];
% 
% title_string = 'Checking that a large error causes fit to not occur';
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [-1 0] + shift_error;
% line_s_start             = 0;
% line_s_end               = 1;
% 
% arc_center_xy            = [0 1];
% arc_radius               = 1;
% arc_vector_start         = [ 0 -1];
% arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% 
% true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% flag_arc_is_first = 0;
% continuity_level = 1;
% [revised_arc1_parameters, revised_arc2_parameters] = fcn_geometry_alignArcToArc(...
%     line_parameters, arc_parameters, flag_arc_is_first, (tolerance), (continuity_level), (fig_num));
% 
% % Check size of results
% assert(isempty(revised_arc1_parameters));
% assert(isempty(revised_arc2_parameters));
% 

%% Basic test 3.1 - an arc nearby the line segment joined with C2 continuity
fig_num = 31;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc1_vector_end           = [ 1  0];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = 0;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [2 1];
arc2_radius               = 3;
arc2_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc2_vector_end           = [ 1  0];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = 0;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

continuity_level = 2;


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc1 is joined to arc2: C2 continuous');


% Check size of results
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_arc1_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));

% Check the arc results
fitted_arc_params = round([revised_arc2_parameters(1,1:3) mod(revised_arc2_parameters(1,4:5),2*pi) revised_arc2_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));

% Check the spiral results
assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_alignArcToArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end