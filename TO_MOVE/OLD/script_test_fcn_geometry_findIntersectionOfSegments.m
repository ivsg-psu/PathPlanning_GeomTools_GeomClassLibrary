% script_test_fcn_geometry_findIntersectionOfSegments
% This is a script to exercise the function: fcn_geometry_findIntersectionOfSegments.m
% This function was written on 2021_06_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Modification history:
%      2021_06_05
%      -- wrote function, adapted from script_test_fcn_geometry_findIntersectionOfSegments.m
%      2024_02_27
%      -- Added test cases for flags 3 and 4. 

clc
close all

%% Simple test 1 - a simple intersection
fprintf(1,'Simple intersection result: \n');
wall_start = [0 10];
wall_end   = [10 10];
sensor_vector_start = [2 1];
sensor_vector_end   = [5 15];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),9.2043));
assert(isequal(round(location,4),[3.9286,10.0000]));


%% Simple test 2 - no intersections
fprintf(1,'No intersection result: \n');
wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type =0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 3 - multiple intersections
fprintf(1,'Multiple intersections result: \n');
wall_start = [0 10; 10 10; 0 6; 10 6];
wall_end = [10 10; 0 6; 10 6; 0 2];

sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 4 - intersection through a vertex
fprintf(1,'Intersection through a vertex result: \n');
wall_start = [0 5; 4 5];
wall_end = [4 5; 8 2];
sensor_vector_start = [4 0];
sensor_vector_end   = [4 8];
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 5 - intersection at start of sensor
fprintf(1,'Intersection at start of sensor result: \n');
path = [0 5; 4 5; 8 2];
wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 5];
sensor_vector_end   = [4 8];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 6 - intersection at end of sensor
fprintf(1,'Intersection at end of sensor result: \n');
path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 5];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 7 - identically overlapping colinear
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
print_more_results(distance,location,path_segments);

%% Simple test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Simple test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Advanced test 1 - intersection beyond a sensor's range with flag
fprintf(1,'Intersection beyond sensor range result: \n');
path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2344;
flag_search_type = 1;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

% Test the negative condition
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];
fig_debugging = 2345;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%   __  __       _ _   _ _    _ _ _
%  |  \/  |     | | | (_) |  | (_) |
%  | \  / |_   _| | |_ _| |__| |_| |_
%  | |\/| | | | | | __| |  __  | | __|
%  | |  | | |_| | | |_| | |  | | | |_
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|
%
%


%% Advanced test 2 - multiple intersections
fprintf(1,'Multiple intersections reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
% sensor_vector_start = [0 0];
% sensor_vector_end   = [5 12];
sensor_vector_start = [10 0];
sensor_vector_end   = [15 12];
fig_debugging = 23488;
flag_search_type = 4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced test 3 - multiple intersections possible, but no hits
fprintf(1,'Multiple intersections possible but no hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 0];
sensor_vector_end   = [0.5 1.2];
fig_debugging = 23499;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced test 4 - multiple intersections possible, but few hits
fprintf(1,'Multiple intersections possible but few hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 0];
sensor_vector_end   = [2.5 6];
fig_debugging = 1010;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _    ____                 _                   _
%  |  \/  |     | | | (_) |  | (_) |  / __ \               | |                 (_)
%  | \  / |_   _| | |_ _| |__| |_| |_| |  | |_   _____ _ __| | __ _ _ __  _ __  _ _ __   __ _
%  | |\/| | | | | | __| |  __  | | __| |  | \ \ / / _ \ '__| |/ _` | '_ \| '_ \| | '_ \ / _` |
%  | |  | | |_| | | |_| | |  | | | |_| |__| |\ V /  __/ |  | | (_| | |_) | |_) | | | | | (_| |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|\____/  \_/ \___|_|  |_|\__,_| .__/| .__/|_|_| |_|\__, |
%                                                                  | |   | |             __/ |
%                                                                  |_|   |_|            |___/


%% Advanced Multihit Overlapping test - identically overlapping colinear
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);






%% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10]; wall_start = path(1:end-1,:); wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _   __  __       _ _ _   _____      _   _
%  |  \/  |     | | | (_) |  | (_) | |  \/  |     | (_) | |  __ \    | | | |
%  | \  / |_   _| | |_ _| |__| |_| |_| \  / |_   _| |_| |_| |__) |_ _| |_| |__
%  | |\/| | | | | | __| |  __  | | __| |\/| | | | | | | __|  ___/ _` | __| '_ \
%  | |  | | |_| | | |_| | |  | | | |_| |  | | |_| | | | |_| |  | (_| | |_| | | |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|_|  |_|\__,_|_|_|\__|_|   \__,_|\__|_| |_|
%
%


%% Advanced Multihit Overlapping test - identically overlapping colinear
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [0 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-2 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [10 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);






%% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [0 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-2 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [-3 10];
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_end = [10 10];
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);



%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [-3 10];
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [13 10];
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Random Multihit
fprintf(1,'random result: \n');

Num_walls = 10;
wall_start = 10*rand(Num_walls,2);
wall_end   = 10*rand(Num_walls,2);
sensor_vector_start = [0 0];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);


%% Advanced Random Single Hit
fprintf(1,'random result: \n');

Num_walls = 10;
wall_start = 10*rand(Num_walls,2);
wall_end   = 10*rand(Num_walls,2);
sensor_vector_start = [0 0];
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Simple test - Extend wall vector (for flag_search_type = 3)
fprintf(1,'No intersection result: \n');
wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type =3;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test - No intersection (for flag_search_type = 3)
fprintf(1,'No intersection result: \n');
wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [4 6];
fig_debugging = 2343;
flag_search_type =3;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test - Extend both the vectors (for flag_search_type = 4)
fprintf(1,'No intersection result: \n');
wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test - Extend both the vectors (for flag_search_type = 4)
fprintf(1,'No intersection result: \n');
wall_start = [-4 10];
wall_end   = [2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [4 6];
fig_debugging = 2343;
flag_search_type =4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced test 1 - intersection beyond a sensor's range with flag (for flag_search_type = 4)
fprintf(1,'Intersection beyond sensor range result: \n');
path = [0 5; 4 5; 8 2];

wall_start = path(1:end-1,:);
wall_end   = path(2:end,:);
sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];
fig_debugging = 2343;
flag_search_type =0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2344;
flag_search_type = 4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

% Test the negative condition
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];
fig_debugging = 2345;
flag_search_type = 0;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2346;
flag_search_type = 4;
[distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end, sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Test of fast implementation mode 

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type, (fig_num));
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
    [distance,location] = ...
    fcn_geometry_findIntersectionOfSegments(...
    wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    flag_search_type, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);



%%
function print_results(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end

%%
function print_more_results(distance,location,path_segments)
fprintf(1,'Distance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end
