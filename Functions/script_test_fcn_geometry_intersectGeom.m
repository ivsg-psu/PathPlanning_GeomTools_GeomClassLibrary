%% script_test_fcn_geometry_intersectGeom
% Exercises the function: fcn_geometry_intersectGeom
% Revision history:
% 2024_05_02 - Aneesh Batchu
% -- wrote the code
% 2024_05_10 - Sean Brennan
% -- added test cases that do not work
% -- added to-do list

% TO-DO
% 2024_05_10 - added by Sean Brennan
% -- make all intersection type figures start with the same number, for
% example all Arc to Arc figures start with 1 (101, 102, 103), all arc to
% line with 2 (201, 202,...), etc.  (this was done for arc-to-arc but needed for all sections)
% -- Add "captions" to each section similar to Arc to Arc so we can find
% bugs easier (this was done for arc-to-arc but needed for all sections)
% -- Fix fail cases where overlaps are not detected in arc-to-arc (search for 191, 192)
%    these are infinite overlap cases. The function should return the FIRST
%    point where the arcs overlap, e.g. the first point on arc2 in this
%    example. There needs to be some type of check for overlapping
%    geometries.
% -- many test cases are missing and are not coded: 
%    arc to circle (0, 1, 2, or infinite intersections)
%    line to line (0, 1, or infinite intersections)
%    line to segment (0, 1, or infinite intersections)
%    segment to line (0, 1, or infinite intersections)
%    segment to segment (0, 1, or infinite intersections)
%    circle to line (etc)
%    circle to segment (etc)

close all

%% check arc to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      _
%     /\              | |            /\
%    /  \   _ __ ___  | |_ ___      /  \   _ __ ___
%   / /\ \ | '__/ __| | __/ _ \    / /\ \ | '__/ __|
%  / ____ \| | | (__  | || (_) |  / ____ \| | | (__
% /_/    \_\_|  \___|  \__\___/  /_/    \_\_|  \___|
%
%
%  _____       _                          _   _
% |_   _|     | |                        | | (_)
%   | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20to%20Arc%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-arc figures start with the number 1

%% Basic Test: Arc to Arc Intersection - no intersections

fig_num = 101;

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

%% Basic Test: Arc to Arc Intersection Case - one intersection

fig_num = 111; 

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

%% Basic Test: Arc to Arc Intersection Case. -- one intersection

fig_num = 112; 

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



%% Basic Test: Arc to Arc Intersection - one intersection

fig_num = 114;

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

%% Basic Test: Arc to Arc Intersection - one intersection - Previous case but the parameters are interchanged


fig_num = 115;

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

%% Basic Test: Arc to Arc Intersection - 2 intersections

fig_num = 116;

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

%% Basic Test: Arc to Arc Intersection - 2 intersections

fig_num = 117;

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

%% Basic Test: Arc to Arc Intersection - 1 intersection at a tangent

fig_num = 118;

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

%% Basic Test: Arc to Arc Intersection - 2 intersections

fig_num = 119;

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

%% Basic Test: Arc to Arc Intersection: 2 intersections - Same case, but the arc parameters are interchanged

fig_num = 120;

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

%% Basic Test: Arc to Arc Intersection - infinite intersections

fig_num = 191;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-70*pi/180) sin(-70*pi/180)];
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

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[0 0]));

%% Basic Test: Arc to Arc Intersection - infinite intersections

fig_num = 192;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-70*pi/180) sin(-70*pi/180)];
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
arc2_parameters = arc1_parameters;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[-3 3]));

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


%% Basic Test: Line to Arc Intersection - BUG (Fixed it)

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
assert(isequal(round(intersection_points,4), [2.9807    2.6602]));

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

%% Basic Test: Line segment to arc Intersection - BUG (Fixed it)

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
assert(isequal(round(intersection_points,4), [1.6583, 0.5000]));

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

%% Basic Test: Arc to Line Intersection - BUG (Fixed it)

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
assert(isequal(round(intersection_points,4), [2.7889 1.8944]));

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
