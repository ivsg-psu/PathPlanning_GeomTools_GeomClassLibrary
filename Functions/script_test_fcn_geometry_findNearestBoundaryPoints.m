% script_test_fcn_geometry_findNearestBoundaryPoints
% Exercises the function: fcn_geometry_findNearestBoundaryPoints
% Revision history:
% 2024_7_25
% Jiabao Zhao wrote the code


%% Test 1, simple example, three points
gridCenters_driven_path = [3 1;2 2;1 3];
gridCenters_non_drivable_grids = [1 1;3 3];
grid_size = 1;
true_boundary_points = [1 1; 3 3];
fig_num = 1;
[true_borders,x_borders,y_borders] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
    gridCenters_non_drivable_grids, gridCenters_driven_path, grid_size,fig_num);
assert(isequal(length(true_borders(:,1)),4));
assert(isequal(length(x_borders),4));
assert(isequal(length(y_borders),4));


%% Test 2  Real Data 
% after running script_test_geometry_updatedSurfaceAnalysis
[true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
    gridCenters_non_drivable_grids, gridCenters_driven_path, grid_size,fig_num);