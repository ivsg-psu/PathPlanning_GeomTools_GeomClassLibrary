%% Introduction to and Purpose of the Geometry Class library
% This is a demonstration script to show the primary functionality of the
% GeomClass library.
%
% This is the explanation of the code that can be found by running
%       script_demo_Geometry.m
% 
% This code repo is typically located at:
%
%   https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu 
% 
% Additional contributers include:
% 2016-2018 - Seth Tau 
% 2023 - Aneesh Batchu
% 2024 - Xinyu Cao


%% Revision History:
% 2023_11_21 - sbrennan@psu.edu
% -- started a demo code set
% 2023_12_27 - sbrennan@psu.edu
% -- switched to using environment settings for checking input parameters
% 2023_12_28 - Aneesh Batchu 
% -- updated Max speed options in all the functions
% -- fast mode updated for all functions
% -- fast mode in fig num updated for all functions
% -- inputs updated for fast mode in all functions
% -- figures not created if fast mode activated
% -- each function has working test scripts to make sure fast mode works 
% -- updated the README.md file.
% 2024_03_14 - sbrennan@psu.edu
% -- updated the Path library dependence to the new function
% 2024_04_11 - sbrennan@psu.edu
% -- add lines and segments to domainBoxByType
% -- used domainBoxByType call in fcn_geometry_fitLinearRegressionFromHoughFit
% -- went through all the scripts and removeed all clc commands
% 2024_04_14 - sbrennan@psu.edu
% -- added joinLineToArc functionality but without spiral joins
% -- added fillArcSequenceTestPoints with true values
% 2024_04_15 - sbrennan@psu.edu
% -- added fcn_geometry_fitSequentialArcs
% -- added methods to compare similar geometries to find the maximum error
% between them (for example, an arc versus an arc or a line versus an arc,
% or a line versus a line), or even points versus points. See functions:
% compareCurves, comparePointsToCurve, and comparePointsToPoints
% -- fixed plotCircle to produce outputs, just like plotArc
% -- ran test scripts to make sure above changes don't break codes
% -- updated the intersection calculation to use the Path
% library function, NOT the geometry library version. Aneesh did this
% earlier but did not add it to the revisions list.
% -- Need to remove long data prints in test scripts
% 2024_04_17 - S. Brennan
% -- added testing times to script_test_all_functions
% -- sped up testing in fcn_geometry_findAngleBetweenAngles, 
% -- sped up testing in fcn_geometry_findPhiConstraints 
% -- added rng seed specification to test script for domainBoxByType since 
%    some results were throwing errors
% -- updated test script and function for
% fcn_geometry_findTangentPointsTwoCircles in prep for region checking
% 2024_04_20 - S. Brennan
% -- added C0, C1, and C2 continuity type join options for alignLineToArc
% 2024_04_25 - Aneesh Batchu
% -- Wrote Assertions to the missing test scripts. 
% -- Also, sped up the test scripts to run faster. 
% 2024_05_06 - Aneesh Batchu
% -- Added codes to calculate intersections between lines and arcs and arcs
% to arcs
% 2024_05_10 - Sean Brennan
% -- Finished arc to arc alignment codes, alignArcToArc
% -- Updated plotGeometry to allow string argument inputs
% 2024_05_14 - Sean Brennan
% -- Finished alignArcSegment
% -- Added flipGeom functionality. Starting to use this in codes.
% 2024_05_16 - Aneesh Batchu
% -- Created "fcn_geometry_fitHoughTransform" to find a cubic polynomial
% hough domain
% 2024_05_17 - Sean Brennan
% -- Fixed fcn_geometry_alignSegmentArc
% 2024_06_14 - Sean Brennan
% -- Fixed fcn_geometry_alignSegmentArc
% 2024_06_19 - 2024_06_21 - Sean Brennan
% -- Reparameterized spiral, line, and segment definitions to be consistent with
% each other and with the other parameter sets
% 2024_06_24 - Aneesh Batchu
% -- created "fcn_geometry_surfaceAnalysis" to classify mapped surfaces
% into drivable and non-drivable surfaces. 
% 2024_07_10 - Sean Brennan
% -- added printing commands to show details of geometric fit, useful for
% debugging geometric fitting steps. See fcn_geometry_printGeometry and
% fcn_geometry_printFitSequences 
% 2024_07_30 - S. Brennan
% -- added spiral approximation function
% -- added function to find nearest boundary points
% -- fixed subtle bugs within fcn_geometry_separatePointsIntoGrids where
% boundary points were not correctly identified in edge cases
% 2024_07_02 - Aneesh Batchu
% -- added a script to find the nearest boundaries of drivable and
% non-drivable surfaces. 

%% To-do items
% 2024_04_15 - S. Brennan
% -- need to check the fcn_geometry_fitSequentialArcs closely. There are
% larger errors at the end points and not sure why.
% 2024_04_11 - Aneesh Batchu
% -- Some functions are missing on the ReadME.md file. It is not up to
% date.
% 2024_04_11 - S. Brennan
% -- Need to fix domain bug in circle regression where the issimplified
% test fails. This is because, in the circle regression, the domain box is
% overlapping itself. This should not be allowed.
% -- create method to remove overlapping domains within Hough fitting (see bug
% above)
% 2025_05_15 - S. Brennan
% -- add spiral and polynomial types to the alignment options


%% Prep the workspace
close all
clc

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools


% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_04_22';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2023_04_22.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'PathClass_v2024_03_14';
library_folders{ith_library} = {'Functions'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/archive/refs/tags/PathClass_v2024_03_14.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'GPSClass_v2023_06_29';
library_folders{ith_library} = {'Functions'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass/archive/refs/tags/GPSClass_v2023_06_29.zip';

% ith_library = ith_library+1;
% library_name{ith_library}    = 'LineFitting_v2023_07_24';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_Association_LineFitting/archive/refs/tags/LineFitting_v2023_07_24.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'FindCircleRadius_v2023_08_02';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_FindCircleRadius/archive/refs/tags/FindCircleRadius_v2023_08_02.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'BreakDataIntoLaps_v2023_08_25';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps/archive/refs/tags/BreakDataIntoLaps_v2023_08_25.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'ParseXODR_v2023_10_23';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_MapTools_ParseXODR/archive/refs/tags/ParseXODR_v2023_10_23.zip';


%% Clear paths and folders, if needed
if 1==0
    clear flag_GeomClass_Folders_Initialized;
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

%% Do we need to set up the work space?
if ~exist('flag_GeomClass_Folders_Initialized','var')
    this_project_folders = {'Functions','Data','LargeData'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);
    flag_GeomClass_Folders_Initialized = 1;
end

%% Set environment flags for input checking
% These are values to set if we want to check inputs or do debugging
% setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','1');
% setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','1');
setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','0');

%% Basic Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Functions
% ____            _        ______                _   _                 
% |  _ \          (_)      |  ____|              | | (_)                
% | |_) | __ _ ___ _  ___  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
% |  _ < / _` / __| |/ __| |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
% | |_) | (_| \__ \ | (__  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
% |____/ \__,_|___/_|\___| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating Euclidean Distance between the points
% Demonstrates fcn_geometry_euclideanPointsToPointsDistance

fig_num = 22;

pt1 = [-1 1 0; 0 0 1; -3 -2 -4];
pt2 = [2 3 4; 4 0 2; -5 3 -2] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);
assert(isequal(round(dist,4), [5.3852; 4.1231; 5.7446]));                                                       

%% fcn_geometry_calcUnitVector - calculates unit vectors in N-D
% Test 1: a basic test
fig_num = 1;
input_vectors = [3 3]; 

unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Test 2: many vectors
fig_num = 2;
input_vectors = randn(10,2); 
unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Test 2: many 3D vectors
fig_num = 3;
input_vectors = randn(10,3); 
unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

%% fcn_geometry_calcOrthogonalVector - calculates orthogonal vectors in N-D coordinates
% Test 1: one vector in 2D
fig_num = 1;
figure(fig_num);
clf;

input_vectors = [3 3]; 

unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, fig_num); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors.*unit_vectors,2);
assert(all(abs(dot_product_sums)<(eps*100)));

% Test 2: many vectors in 2D
fig_num = 2;
figure(fig_num);
clf;

input_vectors = randn(5,2); 

unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, fig_num); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors.*unit_vectors,2);
assert(all(abs(dot_product_sums)<(eps*100)));

% Test 3: a basic test in 3D
fig_num = 3;
figure(fig_num);
clf;

input_vectors = [1/2^0.5 1/2^0.5 0]; 

unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, fig_num); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));

% Test 4: a basic test in 3D, many points
fig_num = 4;
figure(fig_num);
clf;

step = 0.05;
input_vectors = (step:step:1)'.*[3 2 4]; 

unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, fig_num); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));

%% fcn_geometry_separatePointsIntoGrids separates points into grids
% works for 1D, 2D, and 3D points

fig_num = 200;
figure(fig_num);
clf;

inputPoints = 10*rand(300,2);
gridSize = 2;
gridBoundaries = [2 10 2 8]; 

[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),12));
assert(isequal(length(gridCenters(:,1)),12));

%% Circle-related calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Filling%20Test%20Data
 %   _____ _          _                      _       _           _    _____      _            _       _   _                 
 %  / ____(_)        | |                    | |     | |         | |  / ____|    | |          | |     | | (_)                
 % | |     _ _ __ ___| | ___ ______ _ __ ___| | __ _| |_ ___  __| | | |     __ _| | ___ _   _| | __ _| |_ _  ___  _ __  ___ 
 % | |    | | '__/ __| |/ _ \______| '__/ _ \ |/ _` | __/ _ \/ _` | | |    / _` | |/ __| | | | |/ _` | __| |/ _ \| '_ \/ __|
 % | |____| | | | (__| |  __/      | | |  __/ | (_| | ||  __/ (_| | | |___| (_| | | (__| |_| | | (_| | |_| | (_) | | | \__ \
 %  \_____|_|_|  \___|_|\___|      |_|  \___|_|\__,_|\__\___|\__,_|  \_____\__,_|_|\___|\__,_|_|\__,_|\__|_|\___/|_| |_|___/
 % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determining whether an angle lies in b/w two different angles

fig_num = 202;


start_angle_in_radians = 0*pi/180;
change_in_angle = 90*pi/180;
end_angle_in_radians = start_angle_in_radians+change_in_angle;
angles_to_test_in_radians = (start_angle_in_radians:(start_angle_in_radians+2*pi))';
direction = 1;

[isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num));
correct_results = (angles_to_test_in_radians-start_angle_in_radians)<=end_angle_in_radians;
assert(isequal(isAngleBetween,correct_results));

%% Test case for fcn_geometry_findAngleUsing3PointsOnCircle

fig_num = 211;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% Test case for fcn_geometry_findPhiConstraints

fig_num = 121;

p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [-1 1];
[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
fprintf('\t %.2f \t %.2f\n', ...
    mod(phi_start,2*pi)*180/pi, ...
    mod(change,2*pi)*180/pi);

%% Test case for fcn_geometry_findTangentPointFromPointToCircle

fig_num = 2;
centers = [0 0];
radii = 1;
points = [2 3];
cross_prod = -3;
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
    centers,radii,points,cross_prod,fig_num); %#ok<*NASGU>

%% Test case for fcn_geometry_findTangentPointsFromPointToCircle

fig_num = 1;
centers = [0 0];
radii = 1;
points = [2 3];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num); %#ok<*NASGU>

assert(isequal(round(points_tangent,4),[0.9533,-0.3022;-0.6456,0.7637]));

%% Test case for fcn_geometry_findTangentPointsTwoCircles

fig_num = 1551;
centers_start = [0 0];
centers_end   = [0 4];
radii_start   = [0.3]; %#ok<*NBRAK>
radii_end     = [0.2];
flag_inside_or_outside = 0;
[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside, [], [], fig_num) %#ok<*NOPTS,*ASGLU>

assert(isequal(round(points_tangent_start,4),[0.2976,0.0375;-0.2976,0.0375;0.2999,0.0075;-0.2999,0.0075]));
assert(isequal(round(points_tangent_end,4),[-0.1984,3.9750;0.1984,3.9750;0.1999,4.0050;-0.1999,4.005]));

%% Test case for fcn_geometry_findTangentPointTwoCircles

fig_num = 1;
centers_start = [1 1];
centers_end   = [3 1];
radii_start   = [0.5]; %#ok<*NBRAK>
radii_end     = [0.3];
cross_products_start = [ 1];
cross_products_end   = [-1];

[...
    points_tangent_start, ...
    points_tangent_end] ...
    = ...
    fcn_geometry_findTangentPointTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    cross_products_start,...
    cross_products_end,...
    fig_num) %#ok<*NOPTS,*ASGLU>

assert(isequal(round(points_tangent_start,4),[1.2000,0.5417]));
assert(isequal(round(points_tangent_end,4),[2.8800,1.2750]));

%% Test case for fcn_geometry_findIntersectionOfSegments

fprintf(1,'Simple intersection result: \n');
wall_start = [0 10];
wall_end   = [10 10];
sensor_vector_start = [2 1];
sensor_vector_end   = [5 15];
fig_num = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_num);

assert(isequal(round(distance,4),9.2043));
assert(isequal(round(location,4),[3.9286,10.0000]));


%% Calculating a circle's center and radius from 3 points on the circle
fig_num = 102;
points = [0 0; 0.5 4; 1 -1; 4 -3; 6 2; 7 -2; 9 3; 11 3; 15 -0.5];
hold on
figure(1); clf;
for i=1:length(points(:,1))-2
    [centers, radii] = fcn_geometry_circleCenterFrom3Points(points(i:i+2,:),fig_num);
    plot(points(:,1),points(:,2),'r-');
end


%% Plotting a circle
fig_num = 107;
fcn_geometry_plotCircle(centers,radii,'b-',fig_num) 


%% Plotting an arc
% One arc plotting from 0 to 90, positive

fig_num = 108;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
degree_step = []; % Default is to space at 1 degree increments
format = []; % Can be color like 'r.', color matrix like [1 0 1], or complex string

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);

%% calculating which part of a circle is visible
fig_num = 131;
centers = [0 0; 1 4];
radii = [1; 1];
points = [2 3; 3 4];
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
    centers,radii,points,fig_num);

%% Calculating an arc from 2 points and a direction
fig_num = 100;

centers = [0 0; 4 4; 8 10; -6 10];
radii = [1; 2; 4; 3];
start_angles = [90; 0; -90; 45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)]+centers;
end_angles = [45; 135; 180; 0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)]+centers;
cross_products = [-1; 1; -1; 1];

true_angle = start_angles - end_angles;

[angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products,...
    fig_num);

assert(isequal(round(angles,4),round([-pi/4; 3*pi/4; -pi/2; 7*pi/4],4)));

%% Calculating an arc direction from 3 points
fig_num = 102;
points1 = [0 0];
points2 = [-1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,-1))

%% Calculating an angle from points on an arc, from 3 points
% Demos fcn_geometry_arcAngleFrom3Points

fig_num = 103;
Radius = 2;
points = Radius*[[cos(0)    sin(0)]; [cos(pi/4) sin(pi/4)]; [cos(pi/2) sin(pi/2)]];

points1 = points(1,:);
points2 = points(3,:);
points3 = points(2,:);

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,-3*pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,-7*pi/4));
assert(isequal(circle_centers,[0 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,0));

%% Use fcn_geometry_plotArc to generate arc data, in matrix or cell arrays, even without plotting 
% set fig_num to empty after clearing the figure
fig_num = 104;
figure(fig_num); clf;
fig_num = [];

centers = [3 4];
radii = 2; 
start_angle_in_radians = 45 * pi/180;
end_angle_in_radians = 135 * pi/180;
degree_step = []; % Default is 1 degree
format = [];
fig_num = [];

arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format),(fig_num));

% Pull multiple arcs at the same time
centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
degree_step = 5;
format = [];
fig_num = [];

arc_points_cell_array = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format),(fig_num));


fig_num = 104;
figure(fig_num); clf;
figure(fig_num);
hold on; grid on;
axis equal
plot(arc_points_matrix(:,1),arc_points_matrix(:,2),'b.-','MarkerSize',20);
plot(arc_points_cell_array{1}(:,1),arc_points_cell_array{1}(:,2),'k.-','MarkerSize',20);

%% Line Related Calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://patorjk.com/software/taag/#p=display&f=Big&t=Line%20Calculations
%
%  _      _               _____      _            _       _   _
% | |    (_)             / ____|    | |          | |     | | (_)
% | |     _ _ __   ___  | |     __ _| | ___ _   _| | __ _| |_ _  ___  _ __  ___
% | |    | | '_ \ / _ \ | |    / _` | |/ __| | | | |/ _` | __| |/ _ \| '_ \/ __|
% | |____| | | | |  __/ | |___| (_| | | (__| |_| | | (_| | |_| | (_) | | | \__ \
% |______|_|_| |_|\___|  \_____\__,_|_|\___|\__,_|_|\__,_|\__|_|\___/|_| |_|___/
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test case for fcn_geometry_selfCrossProduct

fig_num = 333;
path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);

unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, fig_num); 

%% Test case for fcn_geometry_polarLineFrom2PolarCoords

fig_num = 3333;

start_point_cart = [1 1];
end_point_cart   = [4 3];

[theta1,r1] = cart2pol(start_point_cart(:,1),start_point_cart(:,2));
[theta2,r2]   = cart2pol(end_point_cart(:,1),end_point_cart(:,2));
points = [theta1,r1;theta2,r2];

% Calculate the line
[phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points,fig_num);
axis([0 5 0 5]);
assert(isequal(round(phi,4),-0.9828));
assert(isequal(round(rho,4),-0.2774));

%% Test case for fcn_geometry_flagPointsFurtherFromOriginThanLineSegment

fig_num = 313;
segment_points = [2 3; 4 5];
test_points = [1 4];
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%% Test case for fcn_geometry_flagPointsCloserToOriginThanLineSegment

fig_num = 3113;
segment_points = [2 3; 4 5];
test_points = [1 1.5];
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points,fig_num);
fprintf(1,'Point flags are:\n');
fprintf(1,'\t%.2f\n',point_flags);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Filling%20Test%20Data
%
%  ______ _ _ _ _               _______        _     _____        _
% |  ____(_) | (_)             |__   __|      | |   |  __ \      | |
% | |__   _| | |_ _ __   __ _     | | ___  ___| |_  | |  | | __ _| |_ __ _
% |  __| | | | | | '_ \ / _` |    | |/ _ \/ __| __| | |  | |/ _` | __/ _` |
% | |    | | | | | | | | (_| |    | |  __/\__ \ |_  | |__| | (_| | || (_| |
% |_|    |_|_|_|_|_| |_|\__, |    |_|\___||___/\__| |_____/ \__,_|\__\__,_|
%                        __/ |
%                       |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filling test data for lines and line segments
fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 7 0; 9 5];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

fig_num = 4;
figure(fig_num);
clf;

% It also works in 3D
fig_num = 2;
seed_points = [2 3 0; 4 5 0; 7 0 2; 9 5 3];
M = 10;
sigma = 0.2;

line_3d_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

assert(length(test_points(:,1))>length(seed_points(:,1)))


% Corrupt the results with outliers
fig_num = 3;
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

%% Filling test data for circles
fig_num = 4444;
figure(fig_num);
clf;


circle_center = [3 5];
circle_radius = 2;
M = 10; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num)); %#ok<*NASGU>
axis equal;

%% Filling test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Fill test data for 2 arcs
first_fraction = [0 0.5]; % data from 0 to 50 percent
second_fraction = [0.80 1]; % data from 80 percent to end
N_points = length(onearc_test_points(:,1));

first_fraction_indicies = round(first_fraction*N_points); % find closest indicies
first_fraction_indicies = max([first_fraction_indicies; 1 1],[],1); % Make sure none are below 1
first_fraction_indicies = min([first_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

second_fraction_indicies = round(second_fraction*N_points); % find closest indicies
second_fraction_indicies = max([second_fraction_indicies; 1 1],[],1); % Make sure none are below 1
second_fraction_indicies = min([second_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

twoarc_test_points = ...
    [onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];

corrupted_twoarc_test_points = ...
    [corrupted_onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    corrupted_onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];

% % For debugging
% figure(33838);
% plot(corrupted_twoarc_test_points(:,1),corrupted_twoarc_test_points(:,2),'k.');

%% Filling test data for spheres
fig_num = 6666;
figure(fig_num);
clf;

N_points = 200;
sphere_center = [3 5 0];
sphere_radius = 2;
sigma = 0.02;

sphere_test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num));
assert(length(sphere_test_points(:,1))==N_points);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_sphere_test_points = fcn_geometry_corruptPointsWithOutliers(sphere_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num)); %#ok<*NASGU>
axis equal;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Perfect%20Fits
%
%  _____           __          _     ______ _ _
% |  __ \         / _|        | |   |  ____(_) |
% | |__) |__ _ __| |_ ___  ___| |_  | |__   _| |_ ___
% |  ___/ _ \ '__|  _/ _ \/ __| __| |  __| | | __/ __|
% | |  |  __/ |  | ||  __/ (__| |_  | |    | | |_\__ \
% |_|   \___|_|  |_| \___|\___|\__| |_|    |_|\__|___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perfect line fit
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(seed_points(1:2,:));

%% Perfect circle fit
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

%% Regression fitting of circles
fig_num = 38383;
[best_fit_regression_parameters, best_fit_domain_box, radial_errors, standard_deviation]  = ...
    fcn_geometry_fitCircleRegressionFromHoughFit(circle_test_points(1:3,:),circle_test_points, fig_num); %#ok<*ASGLU>

%% Test case for fcn_geometry_fitVectorToNPoints

fig_num = 4;
A = -3;
B = 2;
C = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*(-A/B) + (-C/B) + 0.2*randn(Npoints,1);
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nFigure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Hough%20Fits
 %  _    _                   _       ______ _ _       
 % | |  | |                 | |     |  ____(_) |      
 % | |__| | ___  _   _  __ _| |__   | |__   _| |_ ___ 
 % |  __  |/ _ \| | | |/ _` | '_ \  |  __| | | __/ __|
 % | |  | | (_) | |_| | (_| | | | | | |    | | |_\__ \
 % |_|  |_|\___/ \__,_|\__, |_| |_| |_|    |_|\__|___/
 %                      __/ |                         
 %                     |___/                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demo Hough line fitting
fig_num = 1111;
figure(fig_num); clf;

transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = 20;

domains = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);



%% Demo Hough line segment fitting
fig_num = 11111;
figure(fig_num); clf;

transverse_tolerance = 0.1;
station_tolerance = 0.4;
points_required_for_agreement = 20;

domains = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);


%% Demo linear pseudo-regression from Hough Line votes
% Show how to do line fitting from Hough votes
fig_num = 111111; % Reuse the previous figure number if you want to overlay
figure(fig_num);clf;

[regression_domain, std_dev_transverse_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(domains{1}, [], fig_num);
fprintf(1,'\n\nFitting results: \n');
fprintf(1,'Expected standard deviation in fit, transverse direction (total least squares), in meters: %.4f\n',sigma);
fprintf(1,'Measured standard deviation in fit, transverse direction (total least squares), in meters: %.4f\n',std_dev_transverse_distance);

%% Demo Hough circle versus arc fitting
% Change the station tolerance to generate different fits

fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% Peform the fits
inputPoints = corrupted_twoarc_test_points;
transverse_tolerance = 0.1;
expected_radii_range = [1 3];
flag_use_permutations = [];

% Use station tolerance low to find only largest arc
station_tolerance = 0.3;
fig_num = 7777;
figure(fig_num); clf;
[best_fitted_Hough_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

% Make station tolerance larger so it finds entire arc, connecting together
% one side to another but not back around thus finding a large arc and not a circle
station_tolerance = 3;
fig_num = 7788;
figure(fig_num); clf;
[best_fitted_Hough_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

% Force fit to a circle by shutting station tolerance off
station_tolerance = [];
fig_num = 7799;
figure(fig_num); clf;
[best_fitted_Hough_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);

%% Demo circle pseudo-regression from Hough Line votes
% Show how to do circle fitting from Hough votes
fig_num = 222222; 
figure(fig_num); clf;

% Grab some data
circle_center = [3 5];
circle_radius = 2;
M = 10; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));

probability_of_corruption = 0.3; % 30 percent of data is bad
magnitude_of_corruption = 3;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

inputPoints = corrupted_circle_test_points;

% Generate Hough votes -  Force fit to a circle by shutting station tolerance off
station_tolerance = [];
flag_use_permutations = [];
transverse_tolerance = 0.1;

[best_fitted_Hough_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);


% Extract out points from the indicies
source_points = inputPoints(best_fit_source_indicies,:);
associated_points_in_domain = inputPoints(best_agreement_indicies,:); 

% Perform the regression fit
[best_fit_regression_parameters, best_fit_domain_box, radial_errors, standard_deviation]  = ...
    fcn_geometry_fitCircleRegressionFromHoughFit(source_points,associated_points_in_domain, fig_num);

% Check the errors
% figure(38383);
% histogram(radial_errors,10)

% The following 3 lines show how to convert the domain box into a
% polyshape, and how to query points using isinterior to identify which
% points are inside the regression domain box
domainPolyShape = polyshape(best_fit_domain_box(:,1),best_fit_domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
IndiciesOfPointsInDomain = isinterior(domainPolyShape,associated_points_in_domain);
best_fit_associated_indicies = find(IndiciesOfPointsInDomain);

fprintf(1,'\n\nResults of circle regression fitting:\n')
fprintf(1,'Actual circle center [X Y] (meters):            %.4f %.4f\n',circle_center(1,1),circle_center(1,2));
fprintf(1,'Hough fitted circle center [X Y] (meters):      %.4f %.4f\n',best_fitted_Hough_parameters(1,1),best_fitted_Hough_parameters(1,2));
fprintf(1,'Regression fitted circle center [X Y] (meters): %.4f %.4f\n',best_fit_regression_parameters(1,1),best_fit_regression_parameters(1,2));
fprintf(1,'Predicted theoretical minimum distance error between actual and fitted center (meters): %.4f\n',sigma/(length(corrupted_circle_test_points(:,1))^0.5));
fprintf(1,'Measured Hough-fit distance error between actual and fitted center (meters):            %.4f\n',sum((circle_center - best_fitted_Hough_parameters(1:2)).^2,2).^0.5);
fprintf(1,'Measured regression-fit distance error between actual and fitted center (meters):       %.4f\n',sum((circle_center - best_fit_regression_parameters(1:2)).^2,2).^0.5);
fprintf(1,'Actual circle radius (meters):            %.4f \n',circle_radius);
fprintf(1,'Fitted Hough circle radius (meters):      %.4f \n',best_fitted_Hough_parameters(1,3));
fprintf(1,'Fitted regression circle radius (meters): %.4f \n',best_fit_regression_parameters(1,3));
fprintf(1,'Radius distance error between actual and fitted (meters) %.4f\n',(circle_radius - best_fit_regression_parameters(1,3)));


%% Demo arc angle measurement
fig_num = 111;
figure(fig_num); clf;

points = twoarc_test_points;
circleCenter = arc_true_circleCenter;
circleRadius = arc_true_circleRadius;
index_source_point = 10;
station_tolerance = 0.5;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);



%% Demo arc pseudo-regression from Hough Line votes
% Show how to do arc fitting from Hough votes
fig_num = 333333; 
figure(fig_num); clf;

% Grab some data
arc_seed_points = [2 3; 4 5; 6 3];

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

inputPoints = corrupted_onearc_test_points;

% Generate Hough votes -  Force fit to a arc by using station tolerance
station_tolerance = 0.5;
flag_use_permutations = [];
transverse_tolerance = 0.1;

[best_fitted_Hough_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, station_tolerance, expected_radii_range, flag_use_permutations, fig_num);


% Extract out points from the indicies
source_points = inputPoints(best_fit_source_indicies,:);
associated_points_in_domain = inputPoints(best_agreement_indicies,:); 

% Perform regression fit

[regression_fit_arc_center_and_radius_and_angles, domain_box, radial_errors, standard_deviation]  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(source_points(1:3,:), associated_points_in_domain, fig_num); 

fprintf(1,'\n\nResults of arc regression fitting:\n')
fprintf(1,'Actual circle center [X Y] (meters): %.4f %.4f\n',arc_true_circleCenter(1,1),arc_true_circleCenter(1,2));
fprintf(1,'Fitted circle center [X Y] (meters): %.4f %.4f\n',regression_fit_arc_center_and_radius_and_angles(1,1),regression_fit_arc_center_and_radius_and_angles(1,2));
fprintf(1,'Predicted max distance error between actual and fitted center (meters)   %.4f\n',sum((arc_true_circleCenter - regression_fit_arc_center_and_radius_and_angles(1:2)).^2,2).^0.5);
fprintf(1,'Measured actual distance error between actual and fitted center (meters) %.4f\n',sigma/(length(circle_test_points(:,1))^0.5));
fprintf(1,'Actual circle radius (meters): %.4f \n',arc_true_circleRadius);
fprintf(1,'Fitted circle radius (meters): %.4f \n',regression_fit_arc_center_and_radius_and_angles(1,3));
fprintf(1,'Radius distance error between actual and fitted (meters) %.4f\n',(arc_true_circleRadius - regression_fit_arc_center_and_radius_and_angles(1,3)));
fprintf(1,'Actual start angle (degrees): %.4f \n',arc_true_start_angle_in_radians*180/pi);
fprintf(1,'Fitted start angle (degrees): %.4f \n',regression_fit_arc_center_and_radius_and_angles(1,4)*180/pi);
fprintf(1,'Actual end angle (degrees):   %.4f \n',arc_true_end_angle_in_radians*180/pi);
fprintf(1,'Fitted end angle (degrees):   %.4f \n',regression_fit_arc_center_and_radius_and_angles(1,5)*180/pi);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Checking%20Fits
%
 %   _____ _               _    _               ______ _ _       
 %  / ____| |             | |  (_)             |  ____(_) |      
 % | |    | |__   ___  ___| | ___ _ __   __ _  | |__   _| |_ ___ 
 % | |    | '_ \ / _ \/ __| |/ / | '_ \ / _` | |  __| | | __/ __|
 % | |____| | | |  __/ (__|   <| | | | | (_| | | |    | | |_\__ \
 %  \_____|_| |_|\___|\___|_|\_\_|_| |_|\__, | |_|    |_|\__|___/
 %                                       __/ |                   
 %                                      |___/                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fcn_geometry_fillEmptyDomainStructure fills an empty domain
emptyDomain = fcn_geometry_fillEmptyDomainStructure;

%% The function fcn_geometry_findAgreementsOfPointsToLineVector checks which points are "hit" by a vector within transverse and station tolerance
rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);


% A basic test of line segment fitting, specifying index-type base_point_index
fig_num = 1;
figure(fig_num); clf;

base_point_index = 1;
base_point = test_points(base_point_index,:);
end_point = [4 5];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point,-1);
transverse_tolerance = 0.2;
station_tolerance = 2;

[agreement_indicies,station_distances] = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

%% fcn_geometry_findArcAgreementIndicies finds which indicies agree with an arc in only station form
seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption)); %, (fig_num));

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];
corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];


fig_num = 234343;
figure(fig_num); clf;

points = inputPoints;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 0.5;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);


%% fcn_geometry_findAgreementsOfPointsToArc finds the indicies that agree with an arc in both radial and station form

fig_num = 233;
figure(fig_num); clf;

inputPoints = corrupted_twoarc_test_points;
base_point_index = 1;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
transverse_tolerance = 0.1;
station_tolerance = 1;
flag_force_circle_fit = [];
threshold_to_check_arc = 1;

[agreement_indicies, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findAgreementsOfPointsToArc(inputPoints, base_point_index, circleCenter, circleRadius, transverse_tolerance, (station_tolerance), (flag_force_circle_fit),(threshold_to_check_arc), (fig_num));

%% Test case for fcn_geometry_findAgreementsOfPointsToCircle 

seed_points = [6 6; 9 3; 6 0];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 8; % Number of points per meter
sigma = 0.02;

outlieronearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

figure(234);
clf;
hold on;
grid on;
grid minor;
axis equal;

corrupted_outlieronearc_test_points= fcn_geometry_corruptPointsWithOutliers(outlieronearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (234));


fig_num = 3333;
figure(fig_num); clf;

points = corrupted_outlieronearc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
transverse_tolerance = 0.2;


agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num));

%% fcn_geometry_plotFitDomains plots the domains of fit
fig_num = 2343;
domains = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement);
fcn_geometry_plotFitDomains(domains, fig_num);

%% Test case for fcn_geometry_findPointsInSequence

fig_num = 4353;

input_distances = [-1 0 3 6 7 8.5 9 10 11.5 13 14 15 16 19 22]';
base_point_index = 6;
station_tolerance = 2;

sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, fig_num);


%% Sequential Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                             _   _       _   ______ _ _   _   _
%  / ____|                           | | (_)     | | |  ____(_) | | | (_)
% | (___   ___  __ _ _   _  ___ _ __ | |_ _  __ _| | | |__   _| |_| |_ _ _ __   __ _
%  \___ \ / _ \/ _` | | | |/ _ \ '_ \| __| |/ _` | | |  __| | | __| __| | '_ \ / _` |
%  ____) |  __/ (_| | |_| |  __/ | | | |_| | (_| | | | |    | | |_| |_| | | | | (_| |
% |_____/ \___|\__, |\__,_|\___|_| |_|\__|_|\__,_|_| |_|    |_|\__|\__|_|_| |_|\__, |
%                 | |                                                           __/ |
%                 |_|                                                          |___/
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Sequential%20Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
% arc_pattern = [...
%     1/20, 15; 
%     0 20;
%     -1/5 10; 
%     0 10;
%     1/15 40; 
%     0 15
%     -1/10 20];

arc_pattern = [...
    1/20, 15; 
    0 20];

M = 10;
sigma = 0.02;

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

% Add noise?
if 1==0
    % Corrupt the results
    probability_of_corruption = 1;
    magnitude_of_corruption = 0.03;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end

% Add outliers?
if 1==0
    % Corrupt the results
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 1;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end


% Initialize the subplots
fig_num_array(1) = fig_num;
fig_num_array(2) = 0;
fig_num_array(3) = 0;
fig_num_array(4) = 0;
figure(fig_num); clf;

% Add starter points (truth) onto subplot 2,1,1
subplot(2,2,1);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the groups of true points
modifiedArcStartIndicies = [trueArcStartIndicies; length(test_points(:,1))];
for ith_plot = 1:length(trueArcStartIndicies(:,1))
    if ~isempty(trueNamedCurveTypes)
        current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,trueNamedCurveTypes{ith_plot},-1);
    else
        current_color = [0 0 0];
    end
    index_range = modifiedArcStartIndicies(ith_plot):modifiedArcStartIndicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end


% Perform the fit forwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));


% Perform the fit backwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 1;
[fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = ...
    fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));


% Compare lengths and parameters
NfitsInSequence = length(fitSequence_points_forward);

% First, make absolutely sure that the number of fits found in the forward
% direction match the same number of fits in the backward direction
if length(fitSequence_points_backward)~=NfitsInSequence
    warning('on','backtrace');
    warning('An error will be thrown at this code location as the fits were directionally different.');
    error('Found different numbers of fits in the fit sequence when comparing forward/backward directions. Code is not able to handle this yet!');
end

for ith_fit = 1:NfitsInSequence
    if ~strcmp(fitSequence_bestFitType_forward{ith_fit},fitSequence_bestFitType_backward{ith_fit})
        warning('on','backtrace');
        warning('An error will be thrown at this code location as the fits were found to be geometrically different.');
        error('Found different geometries as the best fit for the same regions when comparing forward/backward directions. Code is not able to handle this yet!');
    end
end


% Find the probable fit
fitSequence_indicies_matrix_forward = cell2mat(fitSequence_endIndicies_forward)';
fitSequence_indicies_matrix_backward = cell2mat(fitSequence_endIndicies_backward)';
probable_arc_boundary_indicies = round(mean([fitSequence_indicies_matrix_forward fitSequence_indicies_matrix_backward],2));
% probable_arc_boundary_indicies = probable_arc_boundary_indicies(1:end-1,:);

% Print and plot the results
% % Add vertical lines to indicate where the segments are TRUELY changing
% for ith_start = 1:length(arcStartIndicies)
%     plot([arcStartIndicies(ith_start) arcStartIndicies(ith_start)],[-0.1 1.1],'k-','LineWidth',5);
% end
%

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit number:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',ith_fit),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('True start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',trueArcStartIndicies(ith_fit)),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',probable_arc_boundary_indicies(ith_fit)),10));
end
fprintf(1,'\n');

% Add vertical lines to indicate where the segments were identified as
% changing
figure(fig_num_array(1));
subplot(2,2,3);
for ith_start = 1:NfitsInSequence
    
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_start,fitSequence_bestFitType_forward{ith_start},-1);

    plot([probable_arc_boundary_indicies(ith_start) probable_arc_boundary_indicies(ith_start)],[-0.1 1.1],'-','Color',current_color);
end

figure(fig_num_array(1));
subplot(2,2,4);

% Plot the fitted groups of points. If any of the points are mis-labeled,
% there will be one color incorrectly on top of another, for example a red
% point on top of a blue underlying point.
for ith_plot = 1:NfitsInSequence
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,fitSequence_bestFitType_forward{ith_plot},-1);
    index_range = probable_arc_boundary_indicies(ith_plot):probable_arc_boundary_indicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end

%% Connect the fits so that the lines perfectly align with the arcs

fig_num = 23456;
figure(fig_num);clf;

revised_fitSequence_parameters_forward  = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_forward, fitSequence_parameters_forward, fitting_tolerance*2, fig_num);
revised_fitSequence_parameters_backward = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_backward,fitSequence_parameters_backward, fitting_tolerance*2, fig_num);

fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, revised_fitSequence_parameters_forward,(fig_num));
fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, revised_fitSequence_parameters_backward,(fig_num));

subplot(1,2,1);
good_axis_limits = axis;

%% Print the results


NleadCharacters = 20;
threshold = 0.15;

max_forward_error = -inf;
max_backward_error = -inf;
max_averaged_error = -inf;

fprintf(1,'\n\nPARAMETER FIT COMPARISON:\n');
for ith_fit = 1:NfitsInSequence
    fprintf(1,'\n\nFit Sequence Number: %.0d\n', ith_fit); 

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('   '),NleadCharacters));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},1)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('TRUE '),NleadCharacters));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},0)
     
    % fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD'),NleadCharacters));
    % fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit}, fitSequence_parameters_forward{ith_fit},0)
    % 
    % fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE'),NleadCharacters));
    % fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit}, fitSequence_parameters_backward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD REV'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},revised_fitSequence_parameters_forward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE REV'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit},revised_fitSequence_parameters_backward{ith_fit},0)

    averaged_parameters = (revised_fitSequence_parameters_backward{ith_fit} + revised_fitSequence_parameters_forward{ith_fit})/2;

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('AVERAGED'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},averaged_parameters,0)
end
fprintf(1,'\n');

fprintf(1,'Max forward fitting error:  %.3f meters\n',max_forward_error);
fprintf(1,'Max backward fitting error: %.3f meters\n',max_backward_error);
fprintf(1,'Max averaged fitting error: %.3f meters\n',max_averaged_error);


%% Plot the results
comparison_fig_num = 2828;
figure(comparison_fig_num); clf;
hold on;

for ith_fit = 1:NfitsInSequence

    averaged_parameters = (revised_fitSequence_parameters_backward{ith_fit} + revised_fitSequence_parameters_forward{ith_fit})/2;

    sgtitle({'Fit quality',sprintf('Red is %.2fm error, blue is 0m error', threshold)});

    subplot(1,3,1);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_forward{ith_fit}, revised_fitSequence_parameters_forward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_forward_error = max(max_forward_error,max_error);
    title('Forward fitting');
    axis(good_axis_limits)

    subplot(1,3,2);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_backward{ith_fit}, revised_fitSequence_parameters_backward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_backward_error = max(max_backward_error,max_error);
    title('Reverse fitting');
    axis(good_axis_limits)

    subplot(1,3,3);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_backward{ith_fit}, averaged_parameters,...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_averaged_error = max(max_averaged_error,max_error);
    title('Averaged fitting');
    axis(good_axis_limits)


end
%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};
    
    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
%
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
%
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
%
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
%
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add
% % the subfolder path without any sub-subfolder path additions.
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;
    
    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');
        
        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end
        
    end
    
    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);
        
        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end
        
    end
    
    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end
    
    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);
        
        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end
        
        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);
        
        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end
        
        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*'); % BUG FIX - For Macs, must be *, not *.*
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end
        
        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end
        
        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end
        
    end
    
    
    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');
        
        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
    end
    
    
    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.
    
    eval(sprintf('%s = 1;',flag_varname));
end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots
    
    % Nothing to do!
    
    
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies

