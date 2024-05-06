%% script_test_fcn_geometry_intersectGeom
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
circle_radius               = 2;
arc2_vector_start         = [cos(-135*pi/180) sin(-135*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
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

%% Basic Test: Arc to Arc Intersection Case. -- BUG (Fixed it)

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
circle_radius               = 2;
arc2_vector_start         = [cos(-135*pi/180) sin(-135*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(secondFitType,  secondFitType_parameters, firstFitType,  firstFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.0000, 1.2679]));

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
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
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

%% Basic Test: Arc to Arc Intersection - BUG (Fixed it)

fig_num = 14;

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
arc2_center_xy            = [-3 3];
circle_radius               = 2;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.3333    1.1144]));

%% Basic Test: Arc to Arc Intersection - Previous case but the parameters are interchanged


fig_num = 15;

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
arc2_center_xy            = [-3 3];
circle_radius               = 2;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;  
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.3333    1.1144]));

%% Basic Test: Arc to Arc Intersection 

fig_num = 16;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
arc2_center_xy            = [-4 3];
circle_radius               = 2;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(90*pi/180) sin(90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.6250    4.4524]));

%% Basic Test: Arc to Arc Intersection 

fig_num = 17;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
arc2_center_xy            = [-4 3];
circle_radius               = 2;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(90*pi/180) sin(90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.6250    1.5476]));

%% Basic Test: Arc to Arc Intersection 

fig_num = 18;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
arc2_center_xy            = [-1 3];
arc2_radius               = 2;
arc2_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
assert(isequal(round(intersection_points,4),[-3    3]));

%% Basic Test: Arc to Arc Intersection 

fig_num = 19;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
arc2_center_xy            = [-1.5 3];
arc2_radius               = 2;
arc2_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
assert(isequal(round(intersection_points,4),[-2.4167    4.7776]));

%% Basic Test: Arc to Arc Intersection: Same case, but the arc parameters are interchanged

fig_num = 20;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
arc2_center_xy            = [-1.5 3];
arc2_radius               = 2;
arc2_vector_start         = [cos(90*pi/180) sin(90*pi/180)];
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
firstFitType_parameters = arc2_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.4167    4.7776]));

%% Basic Test: Line to Arc Intersection

fig_num = 111;


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
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.2953, 1.0682]));


%% Basic Test: Line to Arc Intersection - BUG

% intersection_points = [NaN NaN]
%
% Reason: the function returns the intersection point that is encountered
% first traversing from the "first" geometry starting at its "start"
% position.
%
% In this case, 

fig_num = 112;


true_line_unit_tangent_vector = [6 2];
true_start_point_xy = [-2 1];
line_s_start             = 0;
line_s_end               = 1;

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;

% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Arc Intersection

fig_num = 113;


true_line_unit_tangent_vector = [1 0];
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
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line to Circle Intersection

fig_num = 221;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 3];

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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[2 2]));
assert(isequal(intersection_points(1,:),[3 3]));
assert(isequal(intersection_points(2,:),[-3 3]));

%% Basic Test: Line to Circle Intersection

fig_num = 222;

true_line_unit_tangent_vector = [0 5];
true_start_point_xy = [-4 3];

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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to arc Intersection

fig_num = 331;

true_line_unit_tangent_vector = [3 0];
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
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to arc Intersection

fig_num = 332;

true_line_unit_tangent_vector = [3 0];
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
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
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

%% Basic Test: Line segment to arc Intersection - BUG

fig_num = 333;

true_line_unit_tangent_vector = [7 0];
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
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
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

%% Basic Test: Line segment to circle Intersection

fig_num = 441;

true_line_unit_tangent_vector = [3 0];
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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to circle Intersection

fig_num = 442;

true_line_unit_tangent_vector = [7 0];
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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[2 2]));
assert(isequal(round(intersection_points(1,:),4), [1.6583, 0.5000]));
assert(isequal(round(intersection_points(2,:),4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to circle Intersection

fig_num = 443;

true_line_unit_tangent_vector = [1 0];
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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
    
firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));


%% Basic Test: Line segment to circle Intersection

fig_num = 444;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];

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
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
    
firstFitType = 'line segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Arc to Line Intersection 

fig_num = 551;

% Fill in arc 1
arc1_center_xy            = [0 3];
circle_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [6 2];
true_start_point_xy = [1 1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.2209, 0.2597]));

%% Basic Test: Arc to Line Intersection - BUG

fig_num = 552;

% Fill in arc 1
arc1_center_xy            = [0 3];
circle_radius               = 3;
arc1_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [4 2];
true_start_point_xy = [1 1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Circle to Arc Intersection Case

fig_num = 11; 

% Fill in arc 1
circle_center_xy            = [-3 3];
circle_radius               = 2;
% arc1_vector_start         = [cos(-150*pi/180) sin(-150*pi/180)];
% arc1_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
% arc1_is_circle            = 0;
% arc1_is_counter_clockwise = 1;
% arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
% arc1_parameters(1,4:5) = arc1_angles;
% arc1_parameters(1,6)   = arc1_is_circle;
% arc1_parameters(1,7)   = arc1_is_counter_clockwise;

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.0000, 1.2679]));


%% Basic Test: Circle to Arc Intersection - BUG?

fig_num = 13;

% Fill in arc 1
circle_center_xy            = [0 3];
circle_radius               = 3;
% arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
% arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
% arc1_is_circle            = 0;
% arc1_is_counter_clockwise = 1;
% arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
% arc1_parameters(1,4:5) = arc1_angles;
% arc1_parameters(1,6)   = arc1_is_circle;
circle_parameters(1,7)   = arc1_is_counter_clockwise;

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));
