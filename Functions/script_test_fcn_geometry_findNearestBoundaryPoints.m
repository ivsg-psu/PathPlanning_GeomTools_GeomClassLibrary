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

fig_num = 1224; 
% after running script_test_geometry_updatedSurfaceAnalysis
[true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
    gridCenters_non_drivable_grids, gridCenters_driven_path, grid_size,fig_num);

%% DATA

% gridCenters_driven_path =
% 
%  -112.3080   55.5933
%  -112.3080   56.3933
%  -111.5080   56.3933
%  -112.3080   57.1933
%  -111.5080   57.1933
%  -111.5080   57.9933
%  -110.7080   57.9933
%  -111.5080   58.7933
%  -110.7080   58.7933
%  -110.7080   59.5933
%  -109.9080   59.5933
%  -110.7080   60.3933
%  -109.9080   60.3933
%  -109.9080   61.1933
%  -109.1080   61.1933
%  -109.1080   61.9933
%  -108.3080   61.9933
%  -109.1080   62.7933
%  -108.3080   62.7933



% gridCenters_non_drivable_grids =
% 
%  -109.1080   53.9933         0
%  -109.1080   54.7933         0
%  -114.7080   55.5933         0
%  -109.1080   55.5933         0
%  -108.3080   55.5933         0
%  -113.9080   56.3933         0
%  -108.3080   56.3933         0
%  -113.9080   57.1933         0
%  -113.1080   57.1933         0
%  -108.3080   57.1933         0
%  -107.5080   57.1933         0
%  -113.1080   57.9933         0
%  -107.5080   57.9933         0
%  -113.1080   58.7933         0
%  -112.3080   58.7933         0
%  -107.5080   58.7933         0
%  -106.7080   58.7933         0
%  -112.3080   59.5933         0
%  -106.7080   59.5933         0
%  -112.3080   60.3933         0
%  -111.5080   60.3933         0
%  -105.9080   60.3933         0
%  -111.5080   61.1933         0
%  -105.9080   61.1933         0
%  -111.5080   61.9933         0
%  -110.7080   61.9933         0
%  -110.7080   62.7933         0
%  -109.9080   63.5933         0

% true_boundary_points =
% 
%  -113.1080   56.7933
%  -112.3080   58.3933
%  -111.5080   59.9933
%  -110.7080   61.5933
%  -109.9080   63.1933
%  -109.1080   55.9933
%  -108.3080   57.5933
%  -107.5080   59.1933
%  -106.7080   59.9933
%  -109.5080   54.7933
%  -109.5080   55.5933
%  -108.7080   56.3933
%  -108.7080   57.1933
%  -107.9080   57.9933
%  -107.9080   58.7933
%  -107.1080   59.5933
%  -106.3080   60.3933
%  -106.3080   61.1933
%  -113.5080   56.3933
%  -112.7080   57.1933
%  -112.7080   57.9933
%  -111.9080   58.7933
%  -111.9080   59.5933
%  -111.1080   60.3933
%  -111.1080   61.1933
%  -110.3080   61.9933
%  -110.3080   62.7933

% 
% grid_size =
% 
%     0.8000

