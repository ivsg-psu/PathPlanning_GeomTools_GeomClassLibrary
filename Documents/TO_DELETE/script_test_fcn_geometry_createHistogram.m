% script_test_fcn_geometry_createHistogram
% Exercises the function: fcn_geometry_createHistogram
% Revision history:
% 2024_7_30
% Jiabao Zhao wrote the code

close all;

%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demonstration case 1: simple example
fig_num = 0001;
figure(fig_num);
clf;

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
fig_num = 0002;
figure(fig_num);
clf;



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
fig_num = 0003;
figure(fig_num);
clf;


% after running Aneesh updata surface anaylysis data
% here the figure number have to be 1 because the the output 
% of the functin is based on the histogram which will break function if the figure number is -1
fig_num = 1;
total_points_in_each_grid_with_points_greater_than_zero = randi([0, 1000], 100, 1);
total_points_in_each_grid_in_the_driven_path = randi([0, 1000], 50, 1);
grid_size = 1;
[point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
    total_points_in_each_grid_with_points_greater_than_zero, grid_size, fig_num);
% no asserts here (real data)