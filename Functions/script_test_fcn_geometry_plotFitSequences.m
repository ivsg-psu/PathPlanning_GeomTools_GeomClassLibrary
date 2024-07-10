%% script_test_fcn_geometry_plotFitSequences
% Exercises the function: fcn_geometry_plotFitSequences
% Revision history:
% 2024_04_12
% -- wrote the code

close all;

%% BASIC test - line plotting
fig_num = 1;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [-1 0];
line_s_start             = 0;
line_s_end               = 1;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3)   = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
line_parameters(1,4)   = line_s_start - line_s_end;


% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_bestFitType{1} = 'line';
fitSequence_bestFitType{2} = 'arc';

fitSequence_parameters{1}  = line_parameters;
fitSequence_parameters{2}  = arc_parameters;


XY_data = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters,[], [], (fig_num));

assert(length(XY_data(1,:))==2);
assert(length(XY_data(:,1))>=0);

% Show the data that was just plotted (to show it's there) by checking the
% figure number
assert(ishandle(fig_num));

% Now plot the data - visually check that it lies atop the results.
plot(XY_data(:,1),XY_data(:,2),'k.');

%% Show how we can pass parameters
fig_num = 2;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [-1 0];
line_s_start             = 0;
line_s_end               = 1;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3)   = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
line_parameters(1,4)   = line_s_start - line_s_end;



% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_bestFitType{1} = 'line';
fitSequence_bestFitType{2} = 'arc';

fitSequence_parameters{1}  = line_parameters;
fitSequence_parameters{2}  = arc_parameters;

segment_length = [];
format_string = sprintf(' ''-'',''Color'',[0.6 0.6 0.6],''LineWidth'',7 ');
XY_data = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters, segment_length, format_string, (fig_num));

assert(length(XY_data(1,:))==2);
assert(length(XY_data(:,1))>=0);

% Show the data that was just plotted (to show it's there) by checking the
% figure number
assert(ishandle(fig_num));

% Now plot the data - visually check that it lies atop the results.
plot(XY_data(:,1),XY_data(:,2),'k.');

%% Show how we can get point data without plotting
fig_num = -1; % NO FIGURE
% figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [-1 0];
line_s_start             = 0;
line_s_end               = 1;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_bestFitType{1} = 'line';
fitSequence_bestFitType{2} = 'arc';

fitSequence_parameters{1}  = line_parameters;
fitSequence_parameters{2}  = arc_parameters;

segment_length = [];
format_string = sprintf(' ''-'',''Color'',[0.6 0.6 0.6],''LineWidth'',7 ');
XY_data = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters, segment_length, format_string, (fig_num));

assert(length(XY_data(1,:))==2);
assert(length(XY_data(:,1))>=0);

% Show the data that was just plotted (to show it's there) by checking the
% figure number
% assert(ishandle(fig_num));

% Now plot the data - visually check that it lies atop the results.
fig_num = 3;
figure(fig_num);
clf;

plot(XY_data(:,1),XY_data(:,2),'k.');
axis equal;

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_plotFitSequences(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end