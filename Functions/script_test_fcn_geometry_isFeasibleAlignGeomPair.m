%% script_test_fcn_geometry_isFeasibleAlignGeomPair
% Tests the function fcn_geometry_isFeasibleAlignGeomPair

% Revision history:
% 2024_05_26 - S Brennan
% -- wrote the code

% The first geometery, if the fig number starts with:
% 1: Circles
% 2: Arcs
% 3: Lines
% 4: Segments

% The second geometery, if the fig number 2nd digit is:
% 1: Circles
% 2: Arcs
% 3: Lines
% 4: Segments

% The third digit is the continuity type

% Fourth number is whether feasible
% 5th and 6th numbers are counts of different tests

%% 2. Arc is the first geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      ______ _          _
%     /\              |  ____(_)        | |
%    /  \   _ __ ___  | |__   _ _ __ ___| |_
%   / /\ \ | '__/ __| |  __| | | '__/ __| __|
%  / ____ \| | | (__  | |    | | |  \__ \ |_
% /_/    \_\_|  \___| |_|    |_|_|  |___/\__|
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Arc%20First
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2.2 Arc is the first geometry, Arc is the second geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%     /\                  /\
%    /  \   _ __ ___     /  \   _ __ ___
%   / /\ \ | '__/ __|   / /\ \ | '__/ __|
%  / ____ \| | | (__   / ____ \| | | (__
% /_/    \_\_|  \___| /_/    \_\_|  \___|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Segment%20Arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C0 Arc to Arc Continuity test - NOT FEASIBLE
fig_num = 220001;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 1;
x_offset = 0.5; 
y_offset = 0; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.5));

%% C0 Segment to Arc Continuity test - FEASIBLE
fig_num = 220101;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 1;
x_offset = 0.1; 
y_offset = 0; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));

%% C1 Arc to Arc Continuity test - NOT FEASIBLE
fig_num = 221001;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 2;
x_offset = 0; 
y_offset = -1; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 1;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),1));

%% C1 Arc to Arc Continuity test - FEASIBLE - CW
fig_num = 221101;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 1;
x_offset = 0.1; 
y_offset = 0; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 1;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));


%% C2 Arc to Arc Continuity test - NOT FEASIBLE
fig_num = 222001;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 2;
x_offset = 0; 
y_offset = 0.1; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 2;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.0707));


%% C2 Arc to Arc Continuity test - FEASIBLE
fig_num = 222101;
figure(fig_num);
clf;



% Fill in arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_radius               = 1;
x_offset = 0; 
y_offset = 0; 
arc1_center_xy            = [x_offset arc1_radius+y_offset];
arc1_vector_start         = [cos(-120*pi/180) sin(-120*pi/180)];
arc1_vector_end           = [cos(- 80*pi/180) sin(- 80*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];

clear arc1_parameters
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

input1_type = 'arc';
input1_parameters = arc1_parameters;

% Fill in arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_radius               = 2;
x_offset = 0; 
y_offset = -0.1; 
arc2_center_xy            = [x_offset arc2_radius+y_offset];
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos(- 20*pi/180) sin(- 20*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];

clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

input2_type = 'arc';
input2_parameters = arc2_parameters;

% Set parameters
continuity_level = 2;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.0707));

%% 2.4 Arc is the first geometry, Segment is the second geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       _____                                 _
%     /\               / ____|                               | |
%    /  \   _ __ ___  | (___   ___  __ _ _ __ ___   ___ _ __ | |_
%   / /\ \ | '__/ __|  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|
%  / ____ \| | | (__   ____) |  __/ (_| | | | | | |  __/ | | | |_
% /_/    \_\_|  \___| |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|
%                                   __/ |
%                                  |___/
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Arc%20Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C0 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 420001;
figure(fig_num);
clf;

% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C0 Segment to Arc Continuity test - FEASIBLE
fig_num = 420101;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));


%% C1 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 421001;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C1 Segment to Arc Continuity test - FEASIBLE
fig_num = 421101;
figure(fig_num);
clf;

% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0.2;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));


%% C2 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 422001;
figure(fig_num);
clf;



% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 2;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C1 Segment to Arc Continuity test - FEASIBLE
fig_num = 422101;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 2;
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0.2;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));



%% 4. Segment is the first geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                                 _     ______ _          _
%  / ____|                               | |   |  ____(_)        | |
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |__   _ _ __ ___| |_
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| |  __| | | '__/ __| __|
%  ____) |  __/ (_| | | | | | |  __/ | | | |_  | |    | | |  \__ \ |_
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__| |_|    |_|_|  |___/\__|
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Segment%20First
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 4.2 Segment is the first geometry, Arc is the second geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                                 _
%  / ____|                               | |       /\
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_     /  \   _ __ ___
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|   / /\ \ | '__/ __|
%  ____) |  __/ (_| | | | | | |  __/ | | | |_   / ____ \| | | (__
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__| /_/    \_\_|  \___|
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Segment%20Arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C0 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 420001;
figure(fig_num);
clf;



% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C0 Segment to Arc Continuity test - FEASIBLE
fig_num = 420101;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 0;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));


%% C1 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 421001;
figure(fig_num);
clf;



% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C1 Segment to Arc Continuity test - FEASIBLE
fig_num = 421101;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 1;
x_offset = 0; 
y_offset = 0.1; 
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0.2;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));


%% C2 Segment to Arc Continuity test - NOT FEASIBLE
fig_num = 422001;
figure(fig_num);
clf;



% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 2;

threshold = 0;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),0));
assert(isequal(round(feasibility_distance,4),0.1));

%% C1 Segment to Arc Continuity test - FEASIBLE
fig_num = 422101;
figure(fig_num);
clf;


% Fill in segment parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-2 0];
segment_length              = 2;

clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


input1_type = 'segment';
input1_parameters = segment_parameters;

% Fill in arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_radius               = 2;
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

input2_type = 'arc';
input2_parameters = arc_parameters;

% Set parameters
continuity_level = 1;

threshold = 0.2;

% Call function
[flag_is_feasible, feasibility_distance] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, (threshold), (fig_num));

% Check sizes
assert(isequal(size(flag_is_feasible),[1 1]));
assert(isequal(size(feasibility_distance),[1 1]));

% Check values
assert(isequal(round(flag_is_feasible,4),1));
assert(isequal(round(feasibility_distance,4),-0.1));



