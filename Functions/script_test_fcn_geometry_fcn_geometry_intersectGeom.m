%% script_test_fcn_geometry_fcn_geometry_intersectGeom
% Exercises the function: fcn_geometry_intersectGeom
% Revision history:
% 2024_05_02 - Aneesh Batchu
% -- wrote the code

close all

%% Basic Test: Arc to Arc Intersection Case

fig_num = 11; 

% Fill in arc 1
arc1_center_xy            = [-3 3];
arc1_radius               = 2;
arc1_vector_start         = [cos(-150*pi/180) sin(-150*pi/180)];
arc1_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
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
arc2_center_xy            = [-1 3];
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

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.0000, 1.2679]));

%% Basic Test: Arc to Arc Intersection Case. 

% In this tes, arc1 and arc2 parameters from the previous test are
% interchaged
% intersection_points = [NaN NaN]
% Reason: the function returns the intersection point that is encountered
% first traversing from the "first" geometry starting at its "start"
% position. 
%
% In this case, the first geometry does not intersect the second geometry.
% The second geometry intersects the first geometry. 

fig_num = 12; 

% Fill in arc 1
arc1_center_xy            = [-3 3];
arc1_radius               = 2;
arc1_vector_start         = [cos(-150*pi/180) sin(-150*pi/180)];
arc1_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
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
arc2_center_xy            = [-1 3];
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

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(secondFitType,  secondFitType_parameters, firstFitType,  firstFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Arc to Arc Intersection 

fig_num = 13;

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

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Arc Intersection

fig_num = 14;


true_line_unit_tangent_vector = [6 2];
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

% Fill in arc 2
arc2_center_xy            = [0 3];
arc2_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
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


firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.2953, 1.0682]));


%% Basic Test: Line to Arc Intersection

fig_num = 15;

shift_error = [0 0.2];

true_line_unit_tangent_vector = [6 2];
true_start_point_xy = [-2 1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

% Fill in arc 2
arc2_center_xy            = [0 3];
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


firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));
