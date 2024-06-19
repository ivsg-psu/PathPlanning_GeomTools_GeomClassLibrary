%% script_test_fcn_geometry_plotGeometry
% Exercises the function: fcn_geometry_plotGeometry
% Revision history:
% 2024_04_12
% -- wrote the code
% 2024_04_15
% -- fixed assertions
% 2024_05_12
% -- added more examples
% 2024_06_19 - Sean Brennan
% -- changed parameter format for spiral to new style:
%            'spiral' - 
%
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
% -- changed parameter format for line to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
% -- changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]


close all;

%% BASIC test - 'none' plotting
fig_num = 10;
figure(fig_num); clf;


segment_length = [];
format_string  = [];
XY_data = fcn_geometry_plotGeometry('none', [],segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC test - circle plotting
fig_num = 1;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_center_xy            = [0 1];
circle_radius               = 1;

circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

segment_length = [];
format_string  = [];

XY_data = fcn_geometry_plotGeometry('circle', circle_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC test - arc plotting
fig_num = 2;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];
format_string  = [];

XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)


%% BASIC test - line plotting
fig_num = 3;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [0 0];

% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_base_point_xy; 
line_parameters(1,3)   = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

segment_length = [];
format_string  = [];
XY_data = fcn_geometry_plotGeometry('line', line_parameters,segment_length, format_string, (fig_num));
axis equal

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)


%% BASIC test - segment plotting
fig_num = 4;
figure(fig_num); clf;


segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-1 0];
segment_length              = 3;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = line_base_point_xy; 
segment_parameters(1,3)   = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


segment_length = [];
format_string  = [];
XY_data = fcn_geometry_plotGeometry('segment', segment_parameters,segment_length, format_string, (fig_num));
axis equal

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC test - spiral plotting
fig_num = 5;
figure(fig_num); clf;
%            'spiral' - 
%
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 

% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
spiral_parameters(1,1:2) = [0 0];
spiral_parameters(1,3)   = pi/4;
spiral_parameters(1,5)   = 0;
spiral_parameters(1,6)   = 20;

segment_length = [];
format_string  = [];
XY_data = fcn_geometry_plotGeometry('spiral', spiral_parameters,segment_length, format_string, (fig_num));
axis equal

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC example  - simple plot string
fig_num = 6;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];
format_string = 'b-.';
XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)


%% BASIC example - complex plot string
fig_num = 7;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];
format_string = sprintf(' ''-'',''Color'',[0.6 0 0],''LineWidth'',7 ');
XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC example - color number string
fig_num = 8;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];
format_string = [1 0 1];
XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC example - linewidth plot string (show it fills color automatically)
fig_num = 7;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];
format_string = sprintf(' ''LineWidth'',7 ');
XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format_string, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_plotGeometry(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end