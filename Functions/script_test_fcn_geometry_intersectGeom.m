%% script_test_fcn_geometry_intersectGeom
% Exercises the function: fcn_geometry_intersectGeom
% Revision history:
% 2024_05_02 - Aneesh Batchu
% - wrote the code
% 2024_05_10 - Sean Brennan
% - added test cases that do not work
% - added to-do list
% 2024_05_15 - Aneesh Batchu
% - modified/created figure numbers based on the format specified below.
% - added test cases for circle to line, circle to segment, segment to
% line and line to segment. 
% 2024_06_19 - Sean Brennan
% - changed parameter format for line to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
% 2024_06_19 - Sean Brennan
% - changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]

% TO-DO
% - Add "captions" to each section similar to Arc to Arc so we can find
% bugs easier (this was done for arc-to-arc but needed for all sections)
% - Fix fail cases where overlaps are not detected in arc-to-arc (search for the fig numbers 13090, 13091)
%    these are infinite overlap cases. The function should return the FIRST
%    point where the arcs overlap, e.g. the first point on arc2 in this
%    example. There needs to be some type of check for overlapping
%    geometries.
% - many test cases are missing and are not coded: 
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
figNum = 11001;
figure(figNum); clf;

% Fill in circle 1
circle1_center_xy            = [-3 3];
circle1_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle1_parameters
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 3];
circle2_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle2_parameters
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Circle Intersection Case - one intersection
figNum = 11101;
figure(figNum); clf;

% Fill in circle 1
circle1_center_xy            = [-3 0];
circle1_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle1_parameters
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 0];
circle2_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle2_parameters
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0.0000, 0.0000]));      


%% Basic Test: Circle to Circle Intersection Case - two intersections
figNum = 11201;
figure(figNum); clf;

% Fill in circle 1
circle1_center_xy            = [-3 0];
circle1_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle1_parameters
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [3 0];
circle2_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle2_parameters
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[2 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    5.1962; 0.0000   -5.1962]));  

%% Basic Test: Circle to Circle Intersection Case - infinite intersections
figNum = 11901;
figure(figNum); clf;

% Fill in circle 1
circle1_center_xy            = [0 0];
circle1_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle1_parameters
circle1_parameters(1,1:2) = circle1_center_xy;
circle1_parameters(1,3)   = circle1_radius;

% Fill in circle 2
circle2_center_xy            = [0 0];
circle2_radius               = 6;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle2_parameters
circle2_parameters(1,1:2) = circle2_center_xy;
circle2_parameters(1,3)   = circle2_radius;

firstFitType = 'circle';
firstFitType_parameters = circle1_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Circle to Arc Intersection Case - no intersections (no-overlapping circles)
figNum = 12001;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Arc Intersection Case - no intersections 
figNum = 12002;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));


%% Basic Test: Circle to Arc Intersection - arc starts on cicle, no intersection? (BUG??)
figNum = 12003;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [0 3];
circle_radius               = 3;
% arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
% arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
% arc1_is_circle            = 0;
% arc1_is_counter_clockwise = 1;
% arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));  

%% Basic Test: Circle to Arc Intersection Case - one intersection
figNum = 12101;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  
                                                                                                                
%% Basic Test: Circle to Arc Intersection Case - two intersections
figNum = 12101;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  

%% Basic Test: Circle to Arc Intersection Case - inf intersections
figNum = 12901;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Line to Circle Intersection - No intersection
figNum = 13001;
figure(figNum); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 5]);
true_start_point_xy = [-4 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Circle Intersection - one intersections
figNum = 13101;
figure(figNum); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [0 6];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0 6]));

%% Basic Test: Line to Circle Intersection - Two intersections
figNum = 13201;
figure(figNum); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Circle to Segment Intersection - No intersection
figNum = 14001;
figure(figNum); clf;

true_line_unit_tangent_vector = [0 5];
true_start_point_xy = [-4 3];

segment_unit_tangent_vector = true_line_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Circle to Segment Intersection - no intersections, segment outside circle
figNum = 14002;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Circle to Segment Intersection - no intersections, segment inside circle

figNum = 14003;
figure(figNum); clf;


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 2;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
    
firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));


%% Basic Test: Circle to Segment Intersection - one intersection, tangent
figNum = 14101;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 6];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0 6]));


%% Basic Test: Circle to Segment Intersection - one intersection, crossing inside to outside
figNum = 14102;
figure(figNum); clf;

true_segment_unit_tangent_vector = [0 1];
true_start_point_xy = [0 3];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 10;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0 6]));

%% Basic Test: Circle to Segment Intersection - Two intersections but outputs the point that is closer to the start point of the segment
figNum = 14103;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 3];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 8;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'circle';
firstFitType_parameters = circle_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-3 3]));


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

close all

%% Basic Test: Circle to Arc Intersection Case - no intersections (no-overlapping circles)
figNum = 21001;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 2;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Circle to Arc Intersection Case - no intersections (overlapping circles)
figNum = 21002;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));


%% Basic Test: Arc to Circle Intersection - arc starts on cicle, no intersection? (BUG??)
figNum = 21003;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [0 3];
circle_radius               = 3;
% arc1_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
% arc1_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
% arc1_is_circle            = 0;
% arc1_is_counter_clockwise = 1;
% arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));  

%% Basic Test: Circle to Arc Intersection Case - one intersection
figNum = 21101;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000]));  
                                                                                                                
%% Basic Test: Circle to Arc Intersection Case - two intersections
figNum = 21101;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [ 0.0000    0.0000])); 

%% Basic Test: Circle to Arc Intersection Case - inf intersections
figNum = 12901;
figure(figNum); clf;

% Fill in arc 1
circle_center_xy            = [-3 0];
circle_radius               = 3;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Arc to Arc Intersection - no intersections
figNum = 22001;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Arc to Arc Intersection Case - one intersection
figNum = 22101; 
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.0000, 1.2679]));

%% Basic Test: Arc to Arc Intersection - one intersection
figNum = 22102;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.3333    1.1144]));

%% Basic Test: Arc to Arc Intersection - 2 intersections, arc1 CW, arc2 CW, intersecting arc2 at 2nd intersection
figNum = 22103;
figure(figNum); clf;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-30*pi/180) sin(-30*pi/180)];
arc1_vector_end           = [cos( 90*pi/180) sin( 90*pi/180)];
arc1_is_circle            = 0;
arc1_is_counter_clockwise = 0;
arc1_angles = [atan2(arc1_vector_start(2),arc1_vector_start(1)); atan2(arc1_vector_end(2),arc1_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_parameters(1,1:2) = arc1_center_xy;
arc1_parameters(1,3)   = arc1_radius;
arc1_parameters(1,4:5) = arc1_angles;
arc1_parameters(1,6)   = arc1_is_circle;
arc1_parameters(1,7)   = arc1_is_counter_clockwise;

% Fill in arc 2
arc2_center_xy            = [-2 0];
circle_radius             = 2;
arc2_vector_start         = [cos(180*pi/180) sin(180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 0;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0    0]));

%% Basic Test: Arc to Arc Intersection - 2 intersections, arc1 CW, arc2 CW, intersecting arc2 at 2nd intersection
figNum = 22104;
figure(figNum); clf;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_end           = [cos(-30*pi/180) sin(-30*pi/180)];
arc1_vector_start         = [cos( 90*pi/180) sin( 90*pi/180)];
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
arc2_center_xy            = [-2 0];
circle_radius             = 2;
arc2_vector_end           = [cos(180*pi/180) sin(180*pi/180)];
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0    0]));



%% Basic Test: Arc to Arc Intersection - 1 intersection at a tangent
figNum = 22105;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-3    3]));

%% Basic Test: Arc to Arc Intersection - 2 intersections
figNum = 22106;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.4167    4.7776]));

%% Basic Test: Arc to Arc Intersection: 2 intersections - Same case, but the arc parameters are interchanged
figNum = 22107;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = arc2_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc2_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc1_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[-2.4167    4.7776]));

%% Basic Test: Arc to Arc Intersection - infinite intersections (OUTPUT: NAN) - Ask Dr. B 
figNum = 22901;
figure(figNum); clf;

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
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);


assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Arc to Arc Intersection - infinite intersections - (BUG??) - Output is NaN

figNum = 22902;
figure(figNum); clf;

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

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points), [1 1]));

%% Check Arc to arc overlapping cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ____                 _                   _                                   _                               _____
%  / __ \               | |                 (_)                 /\              | |            /\               / ____|
% | |  | |_   _____ _ __| | __ _ _ __  _ __  _ _ __   __ _     /  \   _ __ ___  | |_ ___      /  \   _ __ ___  | |     __ _ ___  ___  ___
% | |  | \ \ / / _ \ '__| |/ _` | '_ \| '_ \| | '_ \ / _` |   / /\ \ | '__/ __| | __/ _ \    / /\ \ | '__/ __| | |    / _` / __|/ _ \/ __|
% | |__| |\ V /  __/ |  | | (_| | |_) | |_) | | | | | (_| |  / ____ \| | | (__  | || (_) |  / ____ \| | | (__  | |___| (_| \__ \  __/\__ \
%  \____/  \_/ \___|_|  |_|\__,_| .__/| .__/|_|_| |_|\__, | /_/    \_\_|  \___|  \__\___/  /_/    \_\_|  \___|  \_____\__,_|___/\___||___/
%                               | |   | |             __/ |
%                               |_|   |_|            |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Overlapping%20Arc%20to%20Arc%20Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All arc-to-arc figures start with the number 23

close all

%% Overlap 1: Arc to Arc Intersection - arc1 parameters are exactly same as arc2

figNum = 220001;
figure(figNum); clf;

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
arc2_parameters = arc1_parameters;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points), [1 1]));

%% Overlap 2: Arc to Arc Intersection - arc1 parameters are exactly same as arc2, except arc2_is_counter_clockwise = -1

figNum = 220002;
figure(figNum); clf;

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
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = -1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points), [1 1]));

%% Overlap 3: Arc to Arc Intersection - arc2 is inside arc 1
figNum = 220003;
figure(figNum); clf;

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
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-150*pi/180) sin(-150*pi/180)];
arc2_vector_end           = [cos(-110*pi/180) sin(-110*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points), [1 1]));

%% Overlap 4: Arc to Arc Intersection - arc1 is inside arc 2
figNum = 220004;
figure(figNum); clf;

% Fill in arc 1
arc1_center_xy            = [0 3];
arc1_radius               = 3;
arc1_vector_start         = [cos(-150*pi/180) sin(-150*pi/180)];
arc1_vector_end           = [cos(-110*pi/180) sin(-110*pi/180)];
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
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Arc to Line Intersection - no intersection
figNum = 23001;
figure(figNum); clf;

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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Line Intersection - no intersection
figNum = 23002;
figure(figNum); clf;

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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Line Intersection - one intersection
figNum = 23101;
figure(figNum); clf;

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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.2209, 0.2597]));

%% Basic Test: Arc to Line Intersection  - one intersection
figNum = 23102;

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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [2.7889 1.8944]));

%% Basic Test: Arc to Line Intersection - two intersections - Ask Dr. B
figNum = 23103;
figure(figNum); clf;

% Fill in arc 1
arc1_center_xy            = [3 0];
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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0 0]));

%% Basic Test: Arc to Line Intersection - two intersections
figNum = 23104;
figure(figNum); clf;

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

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);


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

close all

%% Basic Test: Arc to Segment Intersection - zero intersections (outside)
figNum = 24001;
figure(figNum); clf;

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


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [4 0];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - no intersection
figNum = 24002;
figure(figNum); clf;

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



true_segment_unit_tangent_vector = [0 1];
true_start_point_xy = [2 0];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - no intersection
figNum = 24003;
figure(figNum); clf;

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

true_segment_unit_tangent_vector = [-1 1];
true_start_point_xy = [0.1 0.1];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Arc to Segment Intersection - one intersection
figNum = 24101;
figure(figNum); clf;

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


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 0];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;




firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3 0]));

%% Basic Test: Arc to Segment Intersection  - one intersection 
figNum = 24102;

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


true_segment_unit_tangent_vector = [4 2];
true_start_point_xy = [1 1];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 5;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [2.7889 1.8944]));

%% Basic Test: Arc to Segment Intersection  - one intersection  and arc1_is_counter_clockwise = 0
figNum = 24103;

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

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-0.5 -0.1];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


% line_parameters = [1, 0, -0.5, -0.1, 0, 2];

firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-0.4359, -0.1000]));

%% Basic Test: Arc to Segment Intersection - two intersections

figNum = 24104;
figure(figNum); clf;

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


true_segment_unit_tangent_vector = [-1 1];
true_start_point_xy = [-0.5 0];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 5;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.9490, 2.4490]));

%% Basic Test: Arc to Segment Intersection - one intersection
figNum = 24105;
figure(figNum); clf;

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


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0 3]));

%% Basic Test: Arc to Segment Intersection - one intersection EXACTLY at end of arc
figNum = 24106;
figure(figNum); clf;

% Fill in arc 1
arc1_center_xy            = [0 -3];
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

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 0];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



firstFitType = 'arc';
firstFitType_parameters = arc1_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [0 0]));

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

close all

%% Basic Test: Line to Circle Intersection - No intersection
figNum = 31001;
figure(figNum); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 5]);
true_start_point_xy = [-4 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to Circle Intersection - One intersections
figNum = 31101;
figure(figNum); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4),[0 0]));

%% Basic Test: Line to Circle Intersection - Two intersections
figNum = 31201;
figure(figNum); clf;

true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 3];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Line to Arc Intersection - no intersection
figNum = 32001;
figure(figNum); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([1 1]);
true_start_point_xy = [1 8];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 1
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Line to Arc Intersection - no intersection
figNum = 32002;
figure(figNum); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 1]);
true_start_point_xy = [2 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-900*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(all(isnan(intersection_points)));

%% Basic Test: Line to Arc Intersection - One intersection
figNum = 32101;
figure(figNum); clf;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([6 2]);
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));

% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-2.2953, 1.0682]));

%% Basic Test: Line to Arc Intersection - one intersection
figNum = 32102;
figure(figNum); clf;


true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([6 2]);
true_start_point_xy = [-2 1];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));



% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [2.9807    2.6602]));

%% Basic Test: Line to Arc Intersection
figNum = 32103;
figure(figNum); clf;


true_line_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;


firstFitType = 'line';
firstFitType_parameters = line_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line to arc Intersection - two intersections

figNum = 32104;
figure(figNum); clf;

true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-1 1]);
true_start_point_xy = [-0.5 0];

line_unit_tangent_vector = true_line_unit_tangent_vector;
line_base_point_xy       = true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear line_parameters
line_parameters(1,1:2) = line_base_point_xy;
line_parameters(1,3  ) = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));


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

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);


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

close all

%% Basic Test: line to line Intersection - No intersection

figNum = 33001;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));


second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
second_true_start_point_xy = [0 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - line on top of line (WEIRD)
% 
% figNum = 33002;
% figure(figNum); clf;
% 
% first_true_line_unit_tangent_vector = [0 1];
% first_true_start_point_xy = [0 0];
% 
% first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
% first_line_base_point_xy       = first_true_start_point_xy;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% clear first_line_parameters
% first_line_parameters(1,1:2) = first_line_base_point_xy;
% first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));
% 
% 
% 
% second_true_line_unit_tangent_vector = [0 -1];
% second_true_start_point_xy = [0 -2];
% 
% second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
% second_line_base_point_xy       = second_true_start_point_xy;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% clear second_line_parameters
% second_line_parameters(1,1:2) = second_line_base_point_xy;
% second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));
% 
% 
% firstFitType = 'line';
% firstFitType_parameters = first_line_parameters;
% secondFitType = 'line';
% secondFitType_parameters = second_line_parameters;
% 
% intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);
% 
% assert(isequal(size(intersection_points),[1 2]));
% assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - parallel lines, no intersection 
% 
% figNum = 33003;
% figure(figNum); clf;
% 
% first_true_line_unit_tangent_vector = [ 0 1];
% first_true_start_point_xy = [0 0];
% 
% first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
% first_line_base_point_xy       = first_true_start_point_xy;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% clear first_line_parameters
% first_line_parameters(1,1:2) = first_line_base_point_xy;
% first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));
% 
% 
% 
% second_true_line_unit_tangent_vector = [0 -1];
% second_true_start_point_xy = [1 0];
% 
% second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
% second_line_base_point_xy       = second_true_start_point_xy;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% clear second_line_parameters
% second_line_parameters(1,1:2) = second_line_base_point_xy;
% second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));
% 
% 
% firstFitType = 'line';
% firstFitType_parameters = first_line_parameters;
% secondFitType = 'line';
% secondFitType_parameters = second_line_parameters;
% 
% intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);
% 
% assert(isequal(size(intersection_points),[1 2]));
% assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: line to line Intersection - inf intersections due to overlapping lines, but catches first part
figNum = 33004;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));


%% Basic Test: line to line Intersection - one intersection

figNum = 33101;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));

%% Basic Test: Segment to line Intersection - One intersection

figNum = 33102;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;


second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [4 -5];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [4 0]));

%% Basic Test: line to line Intersection - one intersection
figNum = 33103;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: line to line Intersection - infinite intersections (but returns the first intersection point)
figNum = 33104;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Line to Segment Intersection - No intersection

figNum = 34101;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -4];

second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 2;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line to segment Intersection - overlapping - catches first point as one intersection
figNum = 34002;
figure(figNum); clf;

first_true_line_unit_tangent_vector = [1 0];
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = [-1 0];
second_true_start_point_xy = [9 0];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));

%% Basic Test: Line to Segment Intersection - One intersection

figNum = 34101;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]));


%% Basic Test: Line to Segment Intersection - one intersection
figNum = 34102;
figure(figNum); clf;

first_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_line_unit_tangent_vector = first_true_line_unit_tangent_vector;
first_line_base_point_xy       = first_true_start_point_xy;


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_line_parameters
first_line_parameters(1,1:2) = first_line_base_point_xy;
first_line_parameters(1,3  ) = atan2(first_line_unit_tangent_vector(2),first_line_unit_tangent_vector(1));

second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 10;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'line';
firstFitType_parameters = first_line_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: Line to Segment Intersection - infinite intersections (but returns the first intersection point)
figNum = 34103;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 3;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;


second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Line segment to circle Intersection
figNum = 41001;
figure(figNum); clf;


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection
figNum = 41002;
figure(figNum); clf;


true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([1 1]);
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
    
firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection

figNum = 41003;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [0 3];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;
    
firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to circle Intersection - one intersection point
figNum = 41101;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 3;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points(1,:),4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to circle Intersection - Two intersections points
figNum = 41102;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
circle_center_xy            = [0 3];
circle_radius               = 3;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear circle_parameters
circle_parameters(1,1:2) = circle_center_xy;
circle_parameters(1,3)   = circle_radius;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'circle';
secondFitType_parameters = circle_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Line segment to arc Intersection - No intersection
figNum = 42001;
figure(figNum); clf;


true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;


% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Line segment to arc Intersection - One intersection
figNum = 42101;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 3;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [-1.6583, 0.5000]));

%% Basic Test: Line segment to arc Intersection - (Fixed it)
figNum = 42102;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-90*pi/180) sin(-90*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [1.6583, 0.5000]));

%% Basic Test: Line segment to arc Intersection - Two intersections
figNum = 42103;
figure(figNum); clf;

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-4 0.5];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear segment_parameters
segment_parameters(1,1:2) = segment_base_point_xy;
segment_parameters(1,3  ) = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_parameters(1,4)   = segment_length;



% Fill in arc 2
arc2_center_xy            = [0 3];
circle_radius               = 3;
arc2_vector_start         = [cos(-180*pi/180) sin(-180*pi/180)];
arc2_vector_end           = [cos(0*pi/180) sin(0*pi/180)];
arc2_is_circle            = 0;
arc2_is_counter_clockwise = 1;
arc2_angles = [atan2(arc2_vector_start(2),arc2_vector_start(1)); atan2(arc2_vector_end(2),arc2_vector_end(1));];


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear arc2_parameters
arc2_parameters(1,1:2) = arc2_center_xy;
arc2_parameters(1,3)   = circle_radius;
arc2_parameters(1,4:5) = arc2_angles;
arc2_parameters(1,6)   = arc2_is_circle;
arc2_parameters(1,7)   = arc2_is_counter_clockwise;

firstFitType = 'segment';
firstFitType_parameters = segment_parameters;
secondFitType = 'arc';
secondFitType_parameters = arc2_parameters;


intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Segment to line Intersection - No intersection

figNum = 43001;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;


second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to line Intersection - overlapping so infinite points, returns first one
% Segment goes from [0 0] to [0 7], line goes from [9 0] downward
figNum = 43002;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = [1 0];
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;


second_true_line_unit_tangent_vector = [-1 0];
second_true_start_point_xy = [9 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [7 0]));


%% Basic Test: Segment to line Intersection - One intersection

figNum = 43101;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [4 -5];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [4 0]));

%% Basic Test: Segment to line Intersection - one intersection
figNum = 43102;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 9;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]));

%% Basic Test: Segment to line Intersection - infinite intersections (but returns the first intersection point)
figNum = 43103;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 3;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_line_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_line_base_point_xy       = second_true_start_point_xy;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_line_parameters
second_line_parameters(1,1:2) = second_line_base_point_xy;
second_line_parameters(1,3  ) = atan2(second_line_unit_tangent_vector(2),second_line_unit_tangent_vector(1));


firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'line';
secondFitType_parameters = second_line_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

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

close all

%% Basic Test: Segment to segment Intersection - No intersection

figNum = 44001;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 7]);
second_true_start_point_xy = [9 -2];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to segment Intersection - No intersection
figNum = 44002;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 7;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-7 0]);
second_true_start_point_xy = [9 0];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 1;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(isnan(intersection_points),[1 1]));

%% Basic Test: Segment to segment Intersection - one intersection
figNum = 44101;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 9;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;


second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([0 2]);
second_true_start_point_xy = [9 -2];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 4;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [9 0]))

%% Basic Test: Segment to segment Intersection - one intersection
figNum = 44102;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([2 2]);
first_true_start_point_xy = [0 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 9;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([-2 2]);
second_true_start_point_xy = [9 -2];


second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 10;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [3.5, 3.5]))

%% Basic Test: Segment to segment Intersection - infinite intersections (but returns the first intersection point)
figNum = 44103;
figure(figNum); clf;

first_true_segment_unit_tangent_vector = fcn_geometry_calcUnitVector([7 0]);
first_true_start_point_xy = [5 0];

first_segment_base_point_xy       = first_true_start_point_xy;
first_segment_unit_tangent_vector = first_true_segment_unit_tangent_vector;
first_segment_length              = 3;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear first_segment_parameters
first_segment_parameters(1,1:2) = first_segment_base_point_xy;
first_segment_parameters(1,3)   = atan2(first_segment_unit_tangent_vector(2),first_segment_unit_tangent_vector(1));
first_segment_parameters(1,4)   = first_segment_length;



second_true_line_unit_tangent_vector = fcn_geometry_calcUnitVector([5 0]);
second_true_start_point_xy = [0 0];

second_segment_base_point_xy       = second_true_start_point_xy;
second_segment_unit_tangent_vector = second_true_line_unit_tangent_vector;
second_segment_length              = 9;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
clear second_segment_parameters
second_segment_parameters(1,1:2) = second_segment_base_point_xy;
second_segment_parameters(1,3)   = atan2(second_segment_unit_tangent_vector(2),second_segment_unit_tangent_vector(1));
second_segment_parameters(1,4)   = second_segment_length;

firstFitType = 'segment';
firstFitType_parameters = first_segment_parameters;
secondFitType = 'segment';
secondFitType_parameters = second_segment_parameters;

intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, figNum);

assert(isequal(size(intersection_points),[1 2]));
assert(isequal(round(intersection_points,4), [5, 0]))
