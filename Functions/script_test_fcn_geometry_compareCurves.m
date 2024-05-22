%% script_test_fcn_geometry_compareCurves
% Exercises the function: fcn_geometry_compareCurves
% Revision history:
% 2024_04_14
% -- wrote the code
% -- revised from script_test_fcn_geometry_fitVectorToNPoints
%
close all;


%% Test 1: a line as reference curve, another line as the test curve
fig_num = 1;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [0 0];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

reference_curve_type_string = 'line';
reference_curve_parameters  = line_parameters;


line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.1]);
line_base_point_xy       = [0 0];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = [];
curve_test_segment_length = [];

[flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
    (threshold), (curve_test_segment_length), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(points_XY_on_test_curve(1,:))==2);
assert(length(points_XY_on_test_curve(:,1))>=1);
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(points_XY_on_test_curve(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);


%% Test 2: a line as reference curve, another line as the test curve, high threshold
fig_num = 2;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [0 0];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

reference_curve_type_string = 'line';
reference_curve_parameters  = line_parameters;


line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.2]);
line_base_point_xy       = [0 -0.5];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = 2;
curve_test_segment_length = [];

[flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
    (threshold), (curve_test_segment_length), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(points_XY_on_test_curve(1,:))==2);
assert(length(points_XY_on_test_curve(:,1))>=1);
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(points_XY_on_test_curve(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);

%% Test 3: a line as reference curve, another line as the test curve, low threshold
fig_num = 3;
figure(fig_num); clf;

line_unit_tangent_vector = [1 0];
line_base_point_xy       = [0 0];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

reference_curve_type_string = 'line';
reference_curve_parameters  = line_parameters;


line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.2]);
line_base_point_xy       = [0 -0.5];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = 1;
curve_test_segment_length = [];

[flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
    (threshold), (curve_test_segment_length), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(points_XY_on_test_curve(1,:))==2);
assert(length(points_XY_on_test_curve(:,1))>=1);
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(points_XY_on_test_curve(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);


%% Test 201: an arc as reference curve, another line as the test curve, low threshold
fig_num = 201;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [4 -7];
arc_radius               = 10;
arc_vector_start         = ([0 0] - arc_center_xy);
arc_vector_end           = ([8 0] - arc_center_xy);
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

reference_curve_type_string = 'arc';
reference_curve_parameters  = arc_parameters;

line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.2]);
line_base_point_xy       = [0 1];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = 1;
curve_test_segment_length = [];

[flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
    (threshold), (curve_test_segment_length), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(points_XY_on_test_curve(1,:))==2);
assert(length(points_XY_on_test_curve(:,1))>=1);
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(points_XY_on_test_curve(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);



%% Test 202: an arc as reference curve, another line as the test curve, high threshold
fig_num = 202;
figure(fig_num); clf;

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [4 -7];
arc_radius               = 10;
arc_vector_start         = ([0 0] - arc_center_xy);
arc_vector_end           = ([8 0] - arc_center_xy);
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

reference_curve_type_string = 'arc';
reference_curve_parameters  = arc_parameters;

line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.2]);
line_base_point_xy       = [0 1];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = 2;
curve_test_segment_length = [];

[flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
    (threshold), (curve_test_segment_length), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(points_XY_on_test_curve(1,:))==2);
assert(length(points_XY_on_test_curve(:,1))>=1);
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(points_XY_on_test_curve(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);

%% Test of fast implementation mode

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [4 -7];
arc_radius               = 10;
arc_vector_start         = ([0 0] - arc_center_xy);
arc_vector_end           = ([8 0] - arc_center_xy);
arc_is_circle            = 0;
arc_is_counter_clockwise = 0;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6)   = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

reference_curve_type_string = 'arc';
reference_curve_parameters  = arc_parameters;

line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 0.2]);
line_base_point_xy       = [0 1];
line_s_start             = 0;
line_s_end               = 10;


% Fill the line parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

test_curve_type_string = 'line';
test_curve_parameters  = line_parameters;

threshold           = 2;
curve_test_segment_length = [];


% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
        fcn_geometry_compareCurves(...
        reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
        (threshold), (curve_test_segment_length), (fig_num));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
        fcn_geometry_compareCurves(...
        reference_curve_type_string, reference_curve_parameters, test_curve_type_string, test_curve_parameters,...
        (threshold), (curve_test_segment_length), (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_compareCurves:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

assert(averageTimeFast<averageTimeSlow);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_compareCurves(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
