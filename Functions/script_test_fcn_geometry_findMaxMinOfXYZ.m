% script_test_fcn_Geometry_findMaxMinOfXYZ
% Exercises the function: fcn_Geometry_findMaxMinOfXYZ
% Revision history:
% 2024_6_14
% -- wrote the code

close all;


%% Test 1: a basic test
fig_num = 1;
N_points = [1 2 3;2 2 2;1 4 5]

[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_Geometry_findMaxMinOfXYZ(N_points, fig_num); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));
