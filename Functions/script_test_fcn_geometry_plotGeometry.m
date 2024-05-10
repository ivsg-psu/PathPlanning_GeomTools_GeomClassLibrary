%% script_test_fcn_geometry_plotGeometry
% Exercises the function: fcn_geometry_plotGeometry
% Revision history:
% 2024_04_12
% -- wrote the code
% 2024_04_15
% -- fixed assertions

close all;

%% BASIC test - 'none' plotting
fig_num = 10;
figure(fig_num); clf;


segment_length = [];
XY_data = fcn_geometry_plotGeometry('none', [],segment_length, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(fig_num));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% BASIC test - line plotting
fig_num = 1;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [-1 0];
line_s_start             = 0;
line_s_end               = 1;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

segment_length = [];
XY_data = fcn_geometry_plotGeometry('line', line_parameters,segment_length, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(1));
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
arc_parameters(1,7)   = arc_is_counter_clockwise;

segment_length = [];

XY_data = fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, (fig_num));

% Check that a figure opened with this number, and that outputs are right
% sizes
assert(ishandle(2));
assert(length(XY_data(1,:))==2)
assert(length(XY_data(:,1))>1)

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_plotGeometry(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end