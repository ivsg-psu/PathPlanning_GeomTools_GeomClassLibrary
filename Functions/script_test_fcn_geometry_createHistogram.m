% script_test_fcn_geometry_createHistogram
% Exercises the function: fcn_geometry_createHistogram
% Revision history:
% 2024_7_30
% Jiabao Zhao wrote the code
%% Test 1 simple example 
fig_num = 1;
total_points_in_each_grid_with_points_greater_than_zero = [5, 10, 15, 0, 20];
total_points_in_each_grid_in_the_driven_path = [3, 7, 0, 12, 18];
grid_size = 1;
[point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
    total_points_in_each_grid_with_points_greater_than_zero, grid_size, fig_num);
assert(isequal(length(point_density),1));
assert(isequal(length(counts1),20)); % the bin length is 20
assert(isequal(length(counts2),10)); % the bin length is 10
assert(isequal(max(binEdges),20)); % largest number
%% Test 2 simple example
% Example Data
total_points_in_each_grid_with_points_greater_than_zero = randi([0, 1000], 100, 1);
total_points_in_each_grid_in_the_driven_path = randi([0, 1000], 50, 1);
grid_size = 1;
[point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
    total_points_in_each_grid_with_points_greater_than_zero, grid_size, fig_num);
assert(isequal(length(point_density),1));
assert(isequal(length(counts1),20)); % the bin length is 20
assert(isequal(length(counts2),10)); % the bin length is 10
assert(isequal(max(binEdges),1000)); % largest number
%% Test 2 real data 
% after running Aneesh updata surface anaylysis data
% here the figure number have to be 1 because the the output 
% of the functin is baed on the histogram which will break function if the figure number is -1
fig_num = 1;
[point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
    total_points_in_each_grid_with_points_greater_than_zero, grid_size, fig_num);
% no asserts here (real data)