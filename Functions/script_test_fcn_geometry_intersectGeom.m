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
% -- Fix fail cases where overlaps are not detected in arc-to-arc (search for the fig numbers 13090, 13091)
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

% Note:
% figure numbers starting with:
% 1: all circles as the first geometry
% 2: all arcs as the first geometry
% 3: all lines as the first geometry
% 4: all segments as the first geometry
% 5: all ???s as the first geometry

% Note:
% figure numbers with the SECOND number as:
% 1: all circles as the second geometry
% 2: all arcs as the second geometry
% 3: all lines as the second geometry
% 4: all segments as the second geometry
% 5: all ???s as the second geometry

% Third number: The number of expected intersections:
% 0: no intersection cases
% 1: 1 intersection
% 2: 2 intersections
% 9: infinite intersections

% 4th and 5th number: a counter that counts up through the cases in this
% section.

% Example:
% 24206 as a figure number is an arc as the first geometry, intersecting a
% segment as the 2nd geometry, expecting 2 intersections, and the 6th
% example in this situation.

%% check circle to XXXX intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _          _        ______ _          _     _____       _                          _   _
%  / ____(_)        | |      |  ____(_)        | |   |_   _|     | |                        | | (_)
% | |     _ _ __ ___| | ___  | |__   _ _ __ ___| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
% | |    | | '__/ __| |/ _ \ |  __| | | '__/ __| __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
% | |____| | | | (__| |  __/ | |    | | |  \__ \ |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
%  \_____|_|_|  \___|_|\___| |_|    |_|_|  |___/\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circle%20First%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All circle-XXX figures start with the number 1

%% check cirlce to circle intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _          _        _           _____ _          _
%  / ____(_)        | |      | |         / ____(_)        | |
% | |     _ _ __ ___| | ___  | |_ ___   | |     _ _ __ ___| | ___
% | |    | | '__/ __| |/ _ \ | __/ _ \  | |    | | '__/ __| |/ _ \
% | |____| | | | (__| |  __/ | || (_) | | |____| | | | (__| |  __/
%  \_____|_|_|  \___|_|\___|  \__\___/   \_____|_|_|  \___|_|\___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circle%20to%20Circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All circle-to-circle figures start with the number 11

%% Basic Test: Circle to Circle Intersection Case - no intersections
fig_num = 11001;
figure(fig_num); clf;

% Fill in circle 1
circle1_center_xy            = [-3 3];
circle1_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 3];
circle2_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Circle Intersection Case - one intersection
fig_num = 11101;
figure(fig_num); clf;

% Fill in circle 1
circle1_center_xy            = [-3 0];
circle1_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 0];
circle2_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0.0000, 0.0000]));      


%% Basic Test: Circle to Circle Intersection Case - two intersections
fig_num = 11201;
figure(fig_num); clf;

% Fill in circle 1
circle1_center_xy            = [-3 0];
circle1_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 0];
circle2_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[2 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    5.1962; 0.0000   -5.1962]));  

%% Basic Test: Circle to Circle Intersection Case - infinite intersections
fig_num = 11901;
figure(fig_num); clf;

% Fill in circle 1
circle1_center_xy            = [0 0];
circle1_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [0 0];
circle2_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isinf(intersection_points)));

%% check circle to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _          _        _
%  / ____(_)        | |      | |            /\
% | |     _ _ __ ___| | ___  | |_ ___      /  \   _ __ ___
% | |    | | '__/ __| |/ _ \ | __/ _ \    / /\ \ | '__/ __|
% | |____| | | | (__| |  __/ | || (_) |  / ____ \| | | (__
%  \_____|_|_|  \___|_|\___|  \__\___/  /_/    \_\_|  \___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circle%20to%20Arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All circle-to-circle figures start with the number 12

%% Basic Test: Circle to Arc Intersection Case - no intersections (no-overlapping circles)
fig_num = 12001;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [3 0];
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
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Arc Intersection Case - no intersections (overlapping circles)
fig_num = 12002;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos( 0*pi/180) sin( 0*pi/180)];
arc2_vector_end           = [cos(90*pi/180) sin(90*pi/180)];
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
assert(all(isnan(intersection_points)));


%% Basic Test: Circle to Arc Intersection - arc starts on cicle, no intersection? (BUG??)
fig_num = 12003;
figure(fig_num); clf;

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
% circle_parameters(1,7)   = arc1_is_counter_clockwise;

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

%% Basic Test: Circle to Arc Intersection Case - one intersection
fig_num = 12101;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(   0*pi/180) sin(   0*pi/180)];
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
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  
                                                                                                                
%% Basic Test: Circle to Arc Intersection Case - two intersections
fig_num = 12101;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos( 170*pi/180) sin( 170*pi/180)];
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
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  

%% Basic Test: Circle to Arc Intersection Case - inf intersections
fig_num = 12901;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [-3 0];
arc2_radius               = 3;
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos( 170*pi/180) sin( 170*pi/180)];
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
assert(isequal(round(intersection_points,4), [inf, inf]));  
                                                                      
%% check circle to line intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _          _        _          _      _
%  / ____(_)        | |      | |        | |    (_)
% | |     _ _ __ ___| | ___  | |_ ___   | |     _ _ __   ___
% | |    | | '__/ __| |/ _ \ | __/ _ \  | |    | | '_ \ / _ \
% | |____| | | | (__| |  __/ | || (_) | | |____| | | | |  __/
%  \_____|_|_|  \___|_|\___|  \__\___/  |______|_|_| |_|\___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circle%20to%20Line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All circle-to-circle figures start with the number 13

%% Basic Test: Line to Circle Intersection - No intersection
fig_num = 13001;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 5]);
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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Circle Intersection - one intersections
fig_num = 13101;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 6];

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(intersection_points,[0 6]));

%% Basic Test: Line to Circle Intersection - Two intersections
fig_num = 13201;
figure(fig_num); clf;

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[2 2]));
assert(isequal(intersection_points(1,:),[3 3]));
assert(isequal(intersection_points(2,:),[-3 3]));

%% check circle to segment intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _          _        _           _____                                 _
%  / ____(_)        | |      | |         / ____|                               | |
% | |     _ _ __ ___| | ___  | |_ ___   | (___   ___  __ _ _ __ ___   ___ _ __ | |_
% | |    | | '__/ __| |/ _ \ | __/ _ \   \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|
% | |____| | | | (__| |  __/ | || (_) |  ____) |  __/ (_| | | | | | |  __/ | | | |_
%  \_____|_|_|  \___|_|\___|  \__\___/  |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|
%                                                     __/ |
%                                                    |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Circle%20to%20Line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All circle-to-circle figures start with the number 14

%% Basic Test: Circle to Segment Intersection - No intersection
fig_num = 14001;
figure(fig_num); clf;

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Circle to segment to circle Intersection
fig_num = 14002;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([3 0]);
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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection

fig_num = 14003;
figure(fig_num); clf;

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
    
firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));


%% Basic Test: Circle to Segment Intersection - one intersections
fig_num = 14101;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 6];

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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(intersection_points,[0 6]));


%% Basic Test: Circle to Segment Intersection - one intersections
fig_num = 14102;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 6];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = -5;
line_s_end               = 5;


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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(intersection_points,[0 6]));

%% Basic Test: Circle to Segment Intersection - Two intersections but outputs the point that is closer to the start point of the segment
fig_num = 14103;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 8;


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

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(intersection_points,[-3 3]));


%% check arc to XXXX intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      _       ______ _          _     _____       _                          _   _
%     /\              (_)     |  ____(_)        | |   |_   _|     | |                        | | (_)
%    /  \   _ __ ___   _ ___  | |__   _ _ __ ___| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   / /\ \ | '__/ __| | / __| |  __| | | '__/ __| __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  / ____ \| | | (__  | \__ \ | |    | | |  \__ \ |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% /_/    \_\_|  \___| |_|___/ |_|    |_|_|  |___/\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20is%20First%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-XXX figures start with the number 2

%% check arc to circle intersections
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       _           _____ _          _      
%      /\              | |         / ____(_)        | |     
%     /  \   _ __ ___  | |_ ___   | |     _ _ __ ___| | ___ 
%    / /\ \ | '__/ __| | __/ _ \  | |    | | '__/ __| |/ _ \
%   / ____ \| | | (__  | || (_) | | |____| | | | (__| |  __/
%  /_/    \_\_|  \___|  \__\___/   \_____|_|_|  \___|_|\___|
% 
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20to%20Circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-line figures start with the number 21
                                                          
%% Basic Test: Circle to Arc Intersection Case - no intersections (no-overlapping circles)
fig_num = 21001;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [3 0];
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
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Arc Intersection Case - no intersections (overlapping circles)
fig_num = 21002;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos( 0*pi/180) sin( 0*pi/180)];
arc2_vector_end           = [cos(90*pi/180) sin(90*pi/180)];
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
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));


%% Basic Test: Circle to Arc Intersection - arc starts on cicle, no intersection? (BUG??)
fig_num = 21003;
figure(fig_num); clf;

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
% circle_parameters(1,7)   = arc1_is_counter_clockwise;

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
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));  

%% Basic Test: Circle to Arc Intersection Case - one intersection
fig_num = 21101;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(   0*pi/180) sin(   0*pi/180)];
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
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  
                                                                                                                
%% Basic Test: Circle to Arc Intersection Case - two intersections
fig_num = 21101;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [0 2];
arc2_radius               = 2;
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos( 170*pi/180) sin( 170*pi/180)];
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
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000])); 

%% Basic Test: Circle to Arc Intersection Case - inf intersections
fig_num = 12901;
figure(fig_num); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;


% Fill in arc 2
arc2_center_xy            = [-3 0];
arc2_radius               = 3;
arc2_vector_start         = [cos(-100*pi/180) sin(-100*pi/180)];
arc2_vector_end           = [cos( 170*pi/180) sin( 170*pi/180)];
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
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [inf, inf]));  

                                        
%% check arc to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      _
%     /\              | |            /\
%    /  \   _ __ ___  | |_ ___      /  \   _ __ ___
%   / /\ \ | '__/ __| | __/ _ \    / /\ \ | '__/ __|
%  / ____ \| | | (__  | || (_) |  / ____ \| | | (__
% /_/    \_\_|  \___|  \__\___/  /_/    \_\_|  \___|
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20to%20Arc%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-line figures start with the number 22

%% Basic Test: Arc to Arc Intersection - no intersections
fig_num = 22001;
figure(fig_num); clf;

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
fig_num = 22101; 
figure(fig_num); clf;

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

%% Basic Test: Arc to Arc Intersection - one intersection
fig_num = 22102;
figure(fig_num); clf;

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

%% Basic Test: Arc to Arc Intersection - 2 intersections
fig_num = 22103;
figure(fig_num); clf;

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

%% Basic Test: Arc to Arc Intersection - 2 intersections - Same as previous case but the arc parameters are interchanged
fig_num = 22104;
figure(fig_num); clf;

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
fig_num = 22105;
figure(fig_num); clf;

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
fig_num = 22106;
figure(fig_num); clf;

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
fig_num = 22107;
figure(fig_num); clf;

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

%% Basic Test: Arc to Arc Intersection - infinite intersections (OUTPUT: NAN) - Ask Dr. B
fig_num = 22901;
figure(fig_num); clf;

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
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Arc to Arc Intersection - infinite intersections - (BUG??) - Output is NaN

fig_num = 22902;
figure(fig_num); clf;

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
assert(isequal(isnan(intersection_points), [1 1]));

%% check arc to line intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      _          _      _
%     /\              | |        | |    (_)
%    /  \   _ __ ___  | |_ ___   | |     _ _ __   ___
%   / /\ \ | '__/ __| | __/ _ \  | |    | | '_ \ / _ \
%  / ____ \| | | (__  | || (_) | | |____| | | | |  __/
% /_/    \_\_|  \___|  \__\___/  |______|_|_| |_|\___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20to%20Arc%0AIntersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-line figures start with the number 23

%% Basic Test: Arc to Line Intersection - no intersection
fig_num = 23001;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 1]);
true_start_point_xy = [1 8];

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
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Line Intersection - no intersection
fig_num = 23002;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 1]);
true_start_point_xy = [2 0];

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
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Line Intersection - one intersection
fig_num = 23101;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([6 2]);
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

%% Basic Test: Arc to Line Intersection  - one intersection
fig_num = 23102;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([4 2]);
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

%% Basic Test: Arc to Line Intersection - two intersections (BUG??) - Ask Dr. B
fig_num = 23003;
figure(fig_num); clf;

% Fill in arc 1
arc1_center_xy            = [0 3];
circle_radius               = 3;
arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [0 0];

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
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Line Intersection - two intersections
fig_num = 23103;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [-0.5 0];

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
assert(isequal(round(intersection_points,4), [-2.9490, 2.4490]));

%% check arc to segment intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      _           _____                                 _
%     /\              | |         / ____|                               | |
%    /  \   _ __ ___  | |_ ___   | (___   ___  __ _ _ __ ___   ___ _ __ | |_
%   / /\ \ | '__/ __| | __/ _ \   \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|
%  / ____ \| | | (__  | || (_) |  ____) |  __/ (_| | | | | | |  __/ | | | |_
% /_/    \_\_|  \___|  \__\___/  |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|
%                                              __/ |
%                                             |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Arc%20to%20Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-arc figures start with the number 24

%% Basic Test: Arc to Segment Intersection - zero intersections (outside)
fig_num = 24001;
figure(fig_num); clf;

% Fill in arc 1
arc1_center_xy            = [0 0];
circle_radius               = 3;
arc1_vector_start         = [cos( -90*pi/180) sin( -90*pi/180)];
arc1_vector_end           = [cos(  90*pi/180) sin(  90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [4 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - no intersection
fig_num = 24002;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 1]);
true_start_point_xy = [2 0];

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
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - no intersection
fig_num = 24003;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [0.1 0.1];

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
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - one intersection
fig_num = 24101;
figure(fig_num); clf;

% Fill in arc 1
arc1_center_xy            = [0 0];
circle_radius               = 3;
arc1_vector_start         = [cos( -90*pi/180) sin( -90*pi/180)];
arc1_vector_end           = [cos(  90*pi/180) sin(  90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3 0]));

%% Basic Test: Arc to Segment Intersection  - one intersection 
fig_num = 24102;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([4 2]);
true_start_point_xy = [1 1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 5;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [2.7889 1.8944]));

%% Basic Test: Arc to Segment Intersection  - one intersection  and arc1_is_counter_clockwise = 0
fig_num = 24103;

% Fill in arc 1
arc1_center_xy            = [0 -1];
circle_radius               = 1;
% arc1_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
% arc1_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles = [2.9671;1.5708];
% 
% 
% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% arc1_parameters = [0, -1, 1, 2.9671, 1.5708, 0, 0]; 

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-0.5 -0.1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

% line_parameters = [1, 0, -0.5, -0.1, 0, 2];

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-0.4359, -0.1000]));

%% Basic Test: Arc to Segment Intersection - two intersections

fig_num = 24104;
figure(fig_num); clf;

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


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [-0.5 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 5;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.9490, 2.4490]));

%% Basic Test: Arc to Segment Intersection - one intersection
fig_num = 24105;
figure(fig_num); clf;

% Fill in arc 1
arc1_center_xy            = [0 0];
circle_radius               = 3;
arc1_vector_start         = [cos( -90*pi/180) sin( -90*pi/180)];
arc1_vector_end           = [cos( -270*pi/180) sin( -270*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0 3]));

%% Basic Test: Arc to Segment Intersection - one intersection (BUG??) - ask Dr.B
fig_num = 24106;
figure(fig_num); clf;

% Fill in arc 1
arc1_center_xy            = [0 0];
circle_radius               = 3;
arc1_vector_start         = [cos( -90*pi/180) sin( -90*pi/180)];
arc1_vector_end           = [cos( 90*pi/180) sin( 90*pi/180)];  % Check Angle: This is -270 in the previous case
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 1;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;


true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0 3]));

%% check line to XXXX intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _      _              ______ _          _     _____       _                          _   _                 
% | |    (_)            |  ____(_)        | |   |_   _|     | |                        | | (_)                
% | |     _ _ __   ___  | |__   _ _ __ ___| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___ 
% | |    | | '_ \ / _ \ |  __| | | '__/ __| __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
% | |____| | | | |  __/ | |    | | |  \__ \ |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |______|_|_| |_|\___| |_|    |_|_|  |___/\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Line%20First%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All Line-XXXX figures start with the number 3     

%% check line to circle intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _      _              _           _____ _          _
% | |    (_)            | |         / ____(_)        | |
% | |     _ _ __   ___  | |_ ___   | |     _ _ __ ___| | ___
% | |    | | '_ \ / _ \ | __/ _ \  | |    | | '__/ __| |/ _ \
% | |____| | | | |  __/ | || (_) | | |____| | | | (__| |  __/
% |______|_|_| |_|\___|  \__\___/   \_____|_|_|  \___|_|\___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Line%20to%20Circle%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All line-circle figures start with the number 31

%% Basic Test: Line to Circle Intersection - No intersection
fig_num = 31001;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 5]);
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

%% Basic Test: Line to Circle Intersection - One intersections
fig_num = 31101;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0];

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
assert(isequal(intersection_points,[0 0]));

%% Basic Test: Line to Circle Intersection - Two intersections
fig_num = 31201;
figure(fig_num); clf;

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

%% check line to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _      _              _
% | |    (_)            | |            /\
% | |     _ _ __   ___  | |_ ___      /  \   _ __ ___
% | |    | | '_ \ / _ \ | __/ _ \    / /\ \ | '__/ __|
% | |____| | | | |  __/ | || (_) |  / ____ \| | | (__
% |______|_|_| |_|\___|  \__\___/  /_/    \_\_|  \___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Line%20to%20Arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All line-arc figures start with the number 32

%% Basic Test: Line to Arc Intersection - no intersection
fig_num = 32001;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 1]);
true_start_point_xy = [1 8];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

% Fill in arc 1
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
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
assert(all(isnan(intersection_points)));

%% Basic Test: Line to Arc Intersection - no intersection
fig_num = 32002;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 1]);
true_start_point_xy = [2 0];

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
arc2_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
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
assert(all(isnan(intersection_points)));

%% Basic Test: Line to Arc Intersection - One intersection
fig_num = 32101;
figure(fig_num); clf;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([6 2]);
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

%% Basic Test: Line to Arc Intersection - one intersection
fig_num = 32102;
figure(fig_num); clf;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([6 2]);
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
fig_num = 32103;
figure(fig_num); clf;


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

%% Basic Test: Line to arc Intersection - two intersections

fig_num = 32104;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [-0.5 0];

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
arc2_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc2_center_xy;
arc1_parameters(1,3)   = circle_radius;
arc1_parameters(1,4:5) = arc2_angles;
arc1_parameters(1,6)   = arc2_is_circle;
arc1_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);


assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.9490, 2.4490]));


%% check line to line intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _      _              _          _      _
% | |    (_)            | |        | |    (_)
% | |     _ _ __   ___  | |_ ___   | |     _ _ __   ___
% | |    | | '_ \ / _ \ | __/ _ \  | |    | | '_ \ / _ \
% | |____| | | | |  __/ | || (_) | | |____| | | | |  __/
% |______|_|_| |_|\___|  \__\___/  |______|_|_| |_|\___|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Line%20to%20Line%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All line-line figures start with the number 33

%% Basic Test: line to line Intersection - No intersection

fig_num = 33001;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
second_true_start_point_xy = [0 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - No intersection (BUG??) 

fig_num = 33002;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 -7]);
second_true_start_point_xy = [0 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - No intersection 

fig_num = 33003;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 -7]);
second_true_start_point_xy = [1 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - inf intersections (BUG in fcn_Path_findProjectionHitOntoPath??) - Ask Dr.B
fig_num = 33004;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));


%% Basic Test: line to line Intersection - one intersection

fig_num = 33101;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));

%% Basic Test: Segment to line Intersection - One intersection

fig_num = 33102;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [4 -5];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 2;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [4 0]));

%% Basic Test: line to line Intersection - one intersection
fig_num = 33103;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: line to line Intersection - infinite intersections (but returns the first intersection point)
fig_num = 33104;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 3;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [5, 0]));


%% check line to segment intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _      _              _           _____                                 _
% | |    (_)            | |         / ____|                               | |
% | |     _ _ __   ___  | |_ ___   | (___   ___  __ _ _ __ ___   ___ _ __ | |_
% | |    | | '_ \ / _ \ | __/ _ \   \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|
% | |____| | | | |  __/ | || (_) |  ____) |  __/ (_| | | | | | |  __/ | | | |_
% |______|_|_| |_|\___|  \__\___/  |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|
%                                                __/ |
%                                               |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Line%20to%20Line%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All line-segment figures start with the number 34

%% Basic Test: Line to Segment Intersection - No intersection

fig_num = 34101;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -4];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 2;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to segment Intersection - No intersection (BUG in fcn_Path_findProjectionHitOntoPath??) - Ask Dr.B
fig_num = 34002;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Segment Intersection - One intersection

fig_num = 34101;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));


%% Basic Test: Line to Segment Intersection - one intersection
fig_num = 34102;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 10;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: Line to Segment Intersection - infinite intersections (but returns the first intersection point)
fig_num = 34103;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 3;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [5, 0]));


%% check segment to XXXX intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _     ______ _          _     _____       _                          _   _
%  / ____|                               | |   |  ____(_)        | |   |_   _|     | |                        | | (_)
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |__   _ _ __ ___| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| |  __| | | '__/ __| __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  ____) |  __/ (_| | | | | | |  __/ | | | |_  | |    | | |  \__ \ |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__| |_|    |_|_|  |___/\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Segment%20First%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All Segment-XXXX figures start with the number 4

%% check segment to circle intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   _____                                 _     _           _____ _          _      
 %  / ____|                               | |   | |         / ____(_)        | |     
 % | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |_ ___   | |     _ _ __ ___| | ___ 
 %  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| | __/ _ \  | |    | | '__/ __| |/ _ \
 %  ____) |  __/ (_| | | | | | |  __/ | | | |_  | || (_) | | |____| | | | (__| |  __/
 % |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|  \__\___/   \_____|_|_|  \___|_|\___|
 %               __/ |                                                               
 %              |___/                                                                                                                           
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Segment%20to%20Circle%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All segment-circle figures start with the number 41

%% Basic Test: Line segment to circle Intersection
fig_num = 41001;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([3 0]);
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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection
fig_num = 41002;
figure(fig_num); clf;

true_line_unit_tangent_vector = [1 1];
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
    
firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection

fig_num = 41003;
figure(fig_num); clf;

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
    
firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection - one intersection point
fig_num = 41101;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 3;


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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points(1,:),4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to circle Intersection - Two intersections points
fig_num = 41102;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 7;


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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points(1,:),4), [-1.6583, 0.5000]));
% assert(isequal(round(intersection_points(2,:),4), [-1.6583, 0.5000]));

%% check segment to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                                 _     _
%  / ____|                               | |   | |            /\
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |_ ___      /  \   _ __ ___
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| | __/ _ \    / /\ \ | '__/ __|
%  ____) |  __/ (_| | | | | | |  __/ | | | |_  | || (_) |  / ____ \| | | (__
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|  \__\___/  /_/    \_\_|  \___|
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Segment%20to%20Arc%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All segment-arc figures start with the number 42

%% Basic Test: Line segment to arc Intersection - No intersection
fig_num = 42001;
figure(fig_num); clf;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([3 0]);
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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to arc Intersection - One intersection
fig_num = 42101;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([3 0]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 3;


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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to arc Intersection - (Fixed it)
fig_num = 42102;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 7;
% segment_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
% segment_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [1.6583, 0.5000]));

%% Basic Test: Line segment to arc Intersection - Two intersections
fig_num = 42103;
figure(fig_num); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_parameters(1,1:2) = line_unit_tangent_vector;
line_parameters(1,3:4) = line_base_point_xy;
line_parameters(1,5)   = line_s_start;
line_parameters(1,6)   = line_s_end;

% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
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

firstFitType = 'segment';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583    0.5000]));

%% check segment to arc intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                                 _     _          _      _
%  / ____|                               | |   | |        | |    (_)
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |_ ___   | |     _ _ __   ___
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| | __/ _ \  | |    | | '_ \ / _ \
%  ____) |  __/ (_| | | | | | |  __/ | | | |_  | || (_) | | |____| | | | |  __/
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|  \__\___/  |______|_|_| |_|\___|
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Segment%20to%20Line%0A%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All segment-arc figures start with the number 43

%% Basic Test: Segment to line Intersection - No intersection

fig_num = 43001;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to line Intersection - No intersection (BUG in fcn_Path_findProjectionHitOntoPath??) - Ask Dr.B
fig_num = 43002;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));


%% Basic Test: Segment to line Intersection - One intersection

fig_num = 43101;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [4 -5];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 2;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [4 0]));

%% Basic Test: Segment to line Intersection - one intersection
fig_num = 43102;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: Segment to line Intersection - infinite intersections (but returns the first intersection point)
fig_num = 43103;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 3;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [5, 0]));

%% check segment to segment intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                                 _     _           _____                                 _
%  / ____|                               | |   | |         / ____|                               | |
% | (___   ___  __ _ _ __ ___   ___ _ __ | |_  | |_ ___   | (___   ___  __ _ _ __ ___   ___ _ __ | |_
%  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| | __/ _ \   \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __|
%  ____) |  __/ (_| | | | | | |  __/ | | | |_  | || (_) |  ____) |  __/ (_| | | | | | |  __/ | | | |_
% |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|  \__\___/  |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|
%               __/ |                                                   __/ |
%              |___/                                                   |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Segment%20to%20Segment%0A%0A%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All segment-arc figures start with the number 44
%% Basic Test: Segment to segment Intersection - No intersection

fig_num = 44001;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to segment Intersection - No intersection
fig_num = 44002;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 7;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 1;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to segment Intersection - one intersection
fig_num = 44101;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 4;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]))

%% Basic Test: Segment to segment Intersection - one intersection
fig_num = 44102;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 10;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]))

%% Basic Test: Segment to segment Intersection - infinite intersections (but returns the first intersection point)
fig_num = 44103;
figure(fig_num); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;
first_line_s_start             = 0;
first_line_s_end               = 3;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_line_parameters(1,1:2) = first_line_unit_tangent_vector;
first_line_parameters(1,3:4) = first_line_base_point_xy;
first_line_parameters(1,5)   = first_line_s_start;
first_line_parameters(1,6)   = first_line_s_end;

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;
second_line_s_start             = 0;
second_line_s_end               = 9;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_line_parameters(1,1:2) = second_line_unit_tangent_vector;
second_line_parameters(1,3:4) = second_line_base_point_xy;
second_line_parameters(1,5)   = second_line_s_start;
second_line_parameters(1,6)   = second_line_s_end;

firstFitType = 'segment';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [5, 0]))
