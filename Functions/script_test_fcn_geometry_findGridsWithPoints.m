% script_test_fcn_geometry_findGridsWithPoints
% Exercises the function: fcn_geometry_findGridsWithPoints
% Revision history:
% 2024_07_15 - Aneesh Batchu
% -- Wrote the example 
% 2024_07_15 - Jiabao Zhao
% -- Functionlize the code 


%% Test 1 
% figure number
fig_num_LLA = 3001;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% As of now. 
point_density = 60;

[gridIndices_cell_array,total_N_points_in_each_grid,gridCenters, grids_with_zero_points, grids_greater_than_zero_points, gridCenters_zero_point_density, gridCenters_greater_than_zero_point_density] = fcn_geometry_findGridsWithPoints(input_points,grid_size,grid_boundaries,fig_num_LLA);