% script_test_fcn_geometry_classifyGridsintoMappedUnmapped
% Exercises the function: fcn_geometry_classifyGridsintoMappedUnmapped
% Revision history:
% 2024_07_15 - Aneesh Batchu
% -- Wrote the example 
% 2024_07_15 - Jiabao Zhao
% -- Functionlize the code 
%% INPUTS

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

[original_grids_with_required_point_density,gridCenters_low_point_density,gridCenters_required_point_density,current_grids_with_low_point_density,current_grids_with_required_point_density] = fcn_geometry_GridsIntoMappedUnmapped(point_density,total_N_points_in_each_grid,grids_greater_than_zero_points,gridCenters,fig_num_LLA);