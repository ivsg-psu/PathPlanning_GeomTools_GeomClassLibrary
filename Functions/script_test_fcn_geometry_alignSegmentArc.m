%% script_test_fcn_geometry_alignSegmentArc
% Exercises the function: fcn_geometry_alignSegmentArc

% Revision history:
% 2024_04_12 - Sean Brennan
% -- wrote the code
% 2024_04_19 - Sean Brennan
% -- renamed from fcn_geometry_joinLineToArc
% 2024_05_10 - Sean Brennan
% -- added test sections

close all;

% %% check input orientation corrections
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   _____                   _
% %  |_   _|                 | |
% %    | |  _ __  _ __  _   _| |_ ___
% %    | | | '_ \| '_ \| | | | __/ __|
% %   _| |_| | | | |_) | |_| | |_\__ \
% %  |_____|_| |_| .__/ \__,_|\__|___/
% %              | |
% %              |_|
% % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1
% %% Basic test 1.11 - an arc nearby the line segment joined with C1 continuity, arc forwards, line is forwards
% fig_num = 11;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [1 0];
% segment_base_point_xy       = [0.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
% arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
% arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking input corrections, arc forward, line forward: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0    1.0000    1.0000   -2.9671   -1.5708         0    1.0000]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 
% %% Basic test 1.12 - an arc nearby the line segment joined with C1 continuity, arc backwards, line is forwards
% fig_num = 12;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [1 0];
% segment_base_point_xy       = [0.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
% arc_vector_start         = [cos(- 20*pi/180) sin(- 20*pi/180)];
% arc_vector_end           = [cos(-120*pi/180) sin(-120*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking input corrections, arc backward, line forward: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0    1.0000    1.0000   -2.9671   -1.5708         0    1.0000]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 
% %% Basic test 1.13 - an arc nearby the line segment joined with C1 continuity, arc forwards, line is backwards
% fig_num = 13;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [-1 0];
% segment_base_point_xy       = [1.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
% arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
% arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking input corrections, arc forward, line backward: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0    1.0000    1.0000   -2.9671   -1.5708         0    1.0000]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 
% 
% 
% %% Basic test 1.14 - an arc nearby the line segment joined with C1 continuity, arc backward, line is backward
% fig_num = 14;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [-1 0];
% segment_base_point_xy       = [1.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
% arc_vector_start         = [cos(- 20*pi/180) sin(- 20*pi/180)];
% arc_vector_end           = [cos(-120*pi/180) sin(-120*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking input corrections, arc backward, line forward: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0    1.0000    1.0000   -2.9671   -1.5708         0    1.0000]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 
% 
% 
% 
% 
% %% check conversions into St coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   _____ _______                                              _
% %  / ____|__   __|                                            (_)
% % | (___    | |     ______    ___ ___  _ ____   _____ _ __ ___ _  ___  _ __
% %  \___ \   | |    |______|  / __/ _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
% %  ____) |  | |             | (_| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
% % |_____/   |_|              \___\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
% %
% % See: http://patorjk.com/software/taag/#p=display&f=Big&t=ST%20-%20conversion
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The following section checks whether the ST conversion sub-codes are
% % working
% 
% %% Basic test 2.1 - Checking St corrections, arc counter-clockwise
% fig_num = 21;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [1 0];
% segment_base_point_xy       = [0.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
% arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
% arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking St corrections, arc counter-clockwise: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0    1.0000    1.0000   -2.9671   -1.5708         0    1.0000]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 
% %% Basic test 2.2 - Checking St corrections, arc clockwise
% fig_num = 22;
% figure(fig_num); clf;
% 
% tolerance = 0.5; % meters
% 
% % Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% segment_unit_tangent_vector = [1 0];
% segment_base_point_xy       = [0.1 0.2];
% segment_s_start             = 0;
% segment_s_end               = 2;
% 
% segment_parameters(1,1:2) = segment_unit_tangent_vector;
% segment_parameters(1,3:4) = segment_base_point_xy;
% segment_parameters(1,5)   = segment_s_start;
% segment_parameters(1,6)   = segment_s_end;
% 
% 
% 
% % Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_radius               = 1;
% arc_center_xy            = [0 -arc_radius];
% arc_vector_start         = [cos( 170*pi/180) sin( 170*pi/180)];
% arc_vector_end           = [cos(  90*pi/180) sin( 90*pi/180)];
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];
% 
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,6)   = arc_is_circle;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% % Set parameters
% continuity_level = 1;
% 
% % Call function
% [revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
%     segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));
% 
% sgtitle(sprintf('Checking St corrections, arc clockwise: C%.0d continuous',continuity_level));
% 
% 
% % Check sizes
% assert(isequal(size(revised_arc_parameters),[1 7]));
% assert(isequal(size(revised_segment_parameters),[1 6]));
% assert(ischar(revised_intermediate_geometry_join_type));
% assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));
% 
% % Check values
% assert(isequal(round(revised_arc_parameters,4),[0   -1.0000    1.0000    2.9671    1.5708         0         0]));
% assert(isequal(round(revised_segment_parameters,4),[1     0     0     0     0     1]));
% assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
% assert(all(isnan(revised_intermediate_geometry_join_parameters)));
% 

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


%% Basic test 3.11 - checking the (-) cross product, feasible, no intersection
fig_num = 311;
figure(fig_num); clf;

tolerance = 1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = -0.1; 
arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.4359    0.9000    1.0000   -2.0218   -0.3491         0    1.0000]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.12 - checking the (-) cross product, NOT feasible, no intersection
fig_num = 312;
figure(fig_num); clf;

tolerance = 0.1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.13 - checking the (-) cross product, forced NOT feasible, no intersection
fig_num = 313;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.21 - checking the (-) cross product, feasible, intersection, feasible
fig_num = 321;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0.1];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.4000    0.9000    1.0000   -0.9273   -0.3491         0    1.0000]));
assert(isequal(round(revised_segment_parameters,4),[1.0000         0   -2.0000    0.1000         0    3.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.22 - checking the (-) cross product, feasible, intersection, not feasible
fig_num = 322;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0.1];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.31 - checking the (+) cross product, feasible, no intersection
fig_num = 331;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius-0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0; % Arc is clockwise
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000   -1.0000    1.0000    1.5708    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.32 - checking the (+) cross product, NOT feasible, no intersection
fig_num = 332;
figure(fig_num); clf;

tolerance = 0.1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius - 0.3];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.33 - checking the (+) cross product, forced NOT feasible, no intersection
fig_num = 333;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius - 0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 3.41 - checking the (+) cross product, feasible, intersection, feasible
fig_num = 341;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;


% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.5641   -0.9000    1.0000    1.1198    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1.0000    0.0000   -2.0000         0         0    3.0000]));
assert(strcmp(revised_intermediate_geometry_join_type,''));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 3.42 - checking the (+) cross product, feasible, intersection, not feasible
fig_num = 342;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 0;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,''));
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



%% Basic test 4.11 - checking the (-) cross product, feasible, no intersection
fig_num = 411;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000    1.0000    1.0000   -1.5708   -0.3491         0    1.0000]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.12 - checking the (-) cross product, NOT feasible, no intersection
fig_num = 412;
figure(fig_num); clf;

tolerance = 0.1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = +0.1; 
arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.13 - checking the (-) cross product, forced NOT feasible, no intersection
fig_num = 413;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = +0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.21 - checking the (-) cross product, feasible, intersection, feasible
fig_num = 421;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000    1.0000    1.0000   -1.5708   -0.3491         0    1.000]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.22 - checking the (-) cross product, feasible, intersection, not feasible
fig_num = 422;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.31 - checking the (+) cross product, feasible, no intersection
fig_num = 431;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius-0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000   -1.0000    1.0000    1.5708    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.32 - checking the (+) cross product, NOT feasible, no intersection
fig_num = 432;
figure(fig_num); clf;

tolerance = 0.1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius-0.2];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.33 - checking the (+) cross product, forced NOT feasible, no intersection
fig_num = 433;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius-0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 4.41 - checking the (+) cross product, feasible, intersection, feasible
fig_num = 441;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius-0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;


% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000   -1.0000    1.0000    1.5708    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1     0    -2     0     0     2]));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Basic test 4.42 - checking the (+) cross product, feasible, intersection, not feasible
fig_num = 442;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 1;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'segment'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

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


%% Basic test 5.11 - checking the (-) cross product, feasible, no intersection
fig_num = 511;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 1.5;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = 0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[-0.5000    1.1000    1.0000   -0.7877   -0.3491         0    1.0000]));
assert(isequal(round(revised_segment_parameters,4),[1.0000         0   -2.0000         0         0    0.7326]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[1.5662         0   -1.2674    0.0000         0    1.0000]));

%% Basic test 5.12 - checking the (-) cross product, NOT feasible, no intersection
fig_num = 512;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.13 - checking the (-) cross product, forced NOT feasible, no intersection
fig_num = 513;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.21 - checking the (-) cross product, feasible, intersection, feasible
fig_num = 521;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-3 0];
segment_s_start             = 0;
segment_s_end               = 3;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.0000    1.0010    1.0000   -1.4933   -0.3491         0    1.0000]));
assert(isequal(round(revised_segment_parameters,4),[1.0000         0   -3.0000         0         0    2.9225]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[ 0.1550   -0.0000   -0.0775    0.0000         0    1.0000]));

%% Basic test 5.22 - checking the (-) cross product, feasible, intersection, not feasible
fig_num = 522;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-0.5 0.1];

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; y_offset = -0.1; arc_center_xy            = [x_offset arc_radius+y_offset];
arc_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc counter-clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.31 - checking the (+) cross product, feasible, no intersection
fig_num = 531;
figure(fig_num); clf;

tolerance = 0.4; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[  0   -1.0010    1.0000    1.4933    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1.0000         0   -2.0000         0         0    1.9225]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.1550   -0.0000   -0.0775   -0.0000         0   -1.0000]));


%% Basic test 5.32 - checking the (+) cross product, NOT feasible, no intersection
fig_num = 532;
figure(fig_num); clf;

tolerance = 0.1; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.2];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.33 - checking the (+) cross product, forced NOT feasible, no intersection
fig_num = 533;
figure(fig_num); clf;

tolerance = []; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.2];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, forced NOT feasible, no intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));

%% Basic test 5.41 - checking the (+) cross product, feasible, intersection, feasible
fig_num = 541;
figure(fig_num); clf;

tolerance = 0.7; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2.5;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;



% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(isequal(round(revised_arc_parameters,4),[0.5000   -1.0010    1.0000    1.4933    0.1745         0         0]));
assert(isequal(round(revised_segment_parameters,4),[1.0000         0   -2.0000         0         0    2.4225]));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(isequal(round(revised_intermediate_geometry_join_parameters,4),[0.1550   -0.0000    0.4225   -0.0000         0   -1.0000]));


%% Basic test 5.42 - checking the (+) cross product, feasible, intersection, not feasible
fig_num = 542;
figure(fig_num); clf;

tolerance = 0.01; % meters


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_s_start             = 0;
segment_s_end               = 2.5;

segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = segment_s_start;
segment_parameters(1,6)   = segment_s_end;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
arc_center_xy            = [0 -arc_radius+0.1];
arc_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
arc_vector_end           = [cos( 10*pi/180) sin( 10*pi/180)];
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

% Set parameters
continuity_level = 2;

% Call function
[revised_segment_parameters, revised_arc_parameters,  revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
    segment_parameters, arc_parameters, (tolerance), (continuity_level), (fig_num));

sgtitle(sprintf('Checking C%.0d continuous, arc clockwise, NOT feasible, intersection',continuity_level));


% Check sizes
assert(isequal(size(revised_arc_parameters),[1 7]));
assert(isequal(size(revised_segment_parameters),[1 6]));
assert(ischar(revised_intermediate_geometry_join_type));
assert(isequal(size(revised_intermediate_geometry_join_parameters),[1 6]));

% Check values
assert(all(isnan(revised_arc_parameters)));
assert(all(isnan(revised_segment_parameters)));
assert(strcmp(revised_intermediate_geometry_join_type,'spiral'));
assert(all(isnan(revised_intermediate_geometry_join_parameters)));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_alignSegmentArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end