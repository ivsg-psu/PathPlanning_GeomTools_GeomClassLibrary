% script_test_fcn_geometry_findMaxMinOfXYZ
% Exercises the function: fcn_geometry_findMaxMinOfXYZ
% Revision history:
% 2024_6_14
% -- wrote the code

close all;

%% Test 1: 3 by 1
fig_num = 1;
N_points = [1;2;3];
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(N_points, fig_num); 
assert(isequal(Min_x,1))

%% Test: 3 by 2 
fig_num = 1;
N_points = [1 2;2 2;3 4];
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(N_points, fig_num); 
assert(isequal(Min_x,1))


%% Test 3: 3 by 3

fig_num = 1;
N_points = [1 2 3;2 2 2;4 3 2];
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(N_points, fig_num); 
assert(isequal(Min_x,1))

%% Test 4: Random numbers, N by 3

fig_num = 1;
N = 100;
N_points = randn(N,3);
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(N_points, fig_num); 
%% assert(isequal(Min_x,1)) DO NOT COMMENT THIS OUT SINCE IT WILL GIVE A ERROR 