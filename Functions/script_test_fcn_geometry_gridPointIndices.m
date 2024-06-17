%% script_test_fcn_geometry_gridPointIndices
% Exerrcises the function: fcn_geometry_gridPointIndices
% Revision history:
% 2024_6_17
% Jiabao Zhao wrote the code

close all

%% Test 1 
rng(123)
N_points = 1500; 
true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ 9*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data
z = true_z;
points1 = [x, y, z]; 
% fig_num = 1011;
% figure(fig_num)
% 
% view(3)
% plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
% 
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('z')
fig_num = 10333; 
figure(fig_num); clf;
inputPoints = points1;
gridSize = 1;
gridBoundaries = [1 10 1 10 2.5 3.5]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));
output = fcn_geometry_findRepeatedNumber(gridIndices);
cell_array = fcn_geometry_gridPointIndices(gridIndices,fig_num);
assert((length(cell_array{1}))==21);