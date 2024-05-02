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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    1.0000   -2.0944   -1.0311         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3281    2.8683    0.6000    0.7311    1.5708         0    1.0000]));
assert(isequal(round(revised_spiral_join_parameters,4),[1.3217    0.5397    0.5139    2.1421    1.0000    1.6667]));

%% Basic test 5.12 - checking the + to + cross product combination, large to small with shift blocked
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
offset_t = 0;

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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(all(isnan(revised_spiral_join_parameters)));

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
offset_t = 0;

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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.6185         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3455    2.8005    0.6000   -0.3654    1.5708         0    1.0000]));
assert(isequal(round(revised_spiral_join_parameters,4),[0.1898    0.9523    0.8148    2.4202    1.0000    1.6667]));




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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, feasible');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0    3.0000    1.0000   -2.0944   -1.0311         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3281    2.8683    0.6000    0.7311    1.5708         0    1.0000]));
assert(isequal(round(revised_spiral_join_parameters,4),[1.3217    0.5397    0.5139    2.1421    1.0000    1.6667]));

%% Basic test 5.22 - checking the + to + cross product combination, small to large, with shift blocked
fig_num = 522;
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
offset_t = 0;

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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift blocked');

% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc1_parameters)));
assert(all(isnan(revised_arc2_parameters)));
assert(all(isnan(revised_spiral_join_parameters)));

%% Basic test 5.23 - checking the + to + cross product combination, small to large, with shift allowed
fig_num = 523;
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
offset_t = 0;

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


[revised_arc1_parameters, revised_arc2_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignArcToArc(...
    arc1_parameters, arc2_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle('Checking ST conversion: cross-products are + to +, not feasible and shift allowed');


% Check sizes
assert(isequal(size(revised_arc1_parameters),[1 7]));
assert(isequal(size(revised_arc2_parameters),[1 7]));
assert(isequal(size(revised_spiral_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc1_parameters,4),[0.0000    3.0000    1.0000   -2.0944   -0.6185         0    1.0000]));
assert(isequal(round(revised_arc2_parameters,4),[0.3455    2.8005    0.6000   -0.3654    1.5708         0    1.0000]));
assert(isequal(round(revised_spiral_join_parameters,4),[0.1898    0.9523    0.8148    2.4202    1.0000    1.6667]));
















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




%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_alignArcToArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end