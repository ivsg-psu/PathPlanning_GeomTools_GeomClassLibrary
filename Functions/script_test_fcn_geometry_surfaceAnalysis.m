
close all; 

% Revision History
% Aneesh Batchu - 2024_06_21
% -- wrote the code originally

%% Basic Test 1: grid_size = 1. Less data, all the surfaces are mapped: All the surfaces have enough data to fit a plane

% Create data
rng(123)

N_points = 100; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

% Add noise
% true_sigma = 0.1;
% z = true_z + true_sigma*randn(Npoints,1);

points1 = [x, y, z]; 

% Create data

rng(123)

N_points = 100; 

true_parameters = [ 0 0 3]';
points = [2 1].*ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

% Add noise
true_sigma = 0.15;
z = true_z + true_sigma*randn(N_points,1);

points2 = [x, y, z]; 

true_parameters = [ 1 0 1]';
points = [3 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points3 = [x, y, z]; 

fig_num = 12; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points3(:,1),points3(:,2),points3(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 1002; 
figure(fig_num); clf;
input_points = [points1; points2; points3]; 
grid_size = 1;
grid_boundaries = [1 4 1 2 2.5 5.5]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));

XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% % Given x_range and y_range
% x_range = min(gridCenters_mapped_grids(:,1)):grid_size:max(gridCenters_mapped_grids(:,1)); % min_x:gridSize:max_x 
% y_range = min(gridCenters_mapped_grids(:,2)):grid_size:max(gridCenters_mapped_grids(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

x_limits = [min(x_range) max(x_range)];  
y_limits = [0 2]; 
% Calculate boundary points
fig_num = 1002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

% set flag_plot_in_3D = 1 for better understanding of the answers
assert(isequal(drivable_grids, 1))
assert(isequal(non_drivable_grids, [2 6 9]'))
assert(isequal(unmapped_grids, [3 4 5 7 8]'))
assert(isequal(gridCenters_mapped_grids(:,1:2), [1.5, 1.5; 2.5 1.5; 3.5 1.5; 3.5 1.5]))
assert(isequal(drivable_grid_numbers_in_mapped_grids, 1))
assert(isequal(non_drivable_grid_numbers_in_mapped_grids, [2 3 4]'))
assert(isequal(boundary_points, [2, 1.5]))

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

% plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

%% Basic Test 2: grid_size = 0.5. The drivable and drivable surfaces are symmetric. A clear boundary can be seen

rng (123)

N_points = 300;
Ext_Square_Size=3;
Int_Square_Size=1;
fig_num = 22;
diag_flag=1;
noise= 0.2;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);


fig_num = 2002; 
figure(fig_num); clf;
input_points = points; 
grid_size = 0.5;
grid_boundaries = [-3 3 -3 3 -1 1]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));

XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 
x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 2002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%% Basic Test 3: grid size = 1

rng (123)

N_points = 300;
Ext_Square_Size=3;
Int_Square_Size=1;
fig_num = 32;
diag_flag=1;
noise= 0.2;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);

% Surface Analysis

fig_num = 3002; 
figure(fig_num); clf;
input_points = points; 
grid_size = 1;
grid_boundaries = [-3 3 -3 3 -1 1]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));


XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4);

% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 3002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%% Test 4: grid_size = 0.25. The drivable and non-drivable surfaces have no symmetry. The boundary points form a line 

rng (123)

N_points = 2500;
Ext_Square_Size=4;
Int_Square_Size=1;
% ext_point_concentration = 3;
% int_point_concentration = 120; 
fig_num = 42;
diag_flag=1;
noise= 0.2;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);

% points = fcn_geometry_concentricSquaresPointDensity_new(Ext_Square_Size,Int_Square_Size,ext_point_concentration,int_point_concentration,noise,diag_flag,fig_num); 

% % Surface Analysis

fig_num = 4002; 
figure(fig_num); clf;
input_points = points; 
grid_size = 0.4;
% grid_boundaries = [-2 2 -2 2 -1 1]; 
grid_boundaries = [-2 2 -2 2 -1 1];
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1;
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));


XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 4002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%% Test 5: grid_size = 0.25. The drivable and non-drivable surfaces have no symmetry. The boundary points are do not form a line 

% Create data
rng(123)

N_points = 1000; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

% Add noise
% true_sigma = 0.1;
% z = true_z + true_sigma*randn(Npoints,1);

points1 = [x, y, z]; 

% Create data

rng(123)

true_parameters = [ 0 0 3]';
points = [2 1].*ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

% Add noise
true_sigma = 0.15;
z = true_z + true_sigma*randn(N_points,1);

points2 = [x, y, z]; 

true_parameters = [ 1 0 1]';
points = [3 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points3 = [x, y, z]; 

fig_num = 52; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points3(:,1),points3(:,2),points3(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 5002; 
figure(fig_num); clf;
input_points = [points1; points2; points3]; 
grid_size = 0.25;
grid_boundaries = [1 4 1 2 2.5 5.5]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));


XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 5002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%% Test 6:

% Create data

rng(123)

N_points = 300; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ [4,1].*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

% Add noise
% true_sigma = 0.1;
% z = true_z + true_sigma*randn(Npoints,1);

points1 = [x, y, z]; 

% Create data

rng(123)

N_points = 300; 

true_parameters = [ 0 0 3]';
points = [5 1].*ones(N_points,1)+ [4,1].*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

% Add noise
true_sigma = 0.2;
z = true_z + true_sigma*randn(N_points,1);

points2 = [x, y, z]; 

% true_parameters = [ 1 0 1]';
% points = [3 1].*ones(N_points,2) + rand(N_points,2);
% x = points(:,1);
% y = points(:,2);
% true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data
% 
% z = true_z;
% 
% points3 = [x, y, z]; 

fig_num = 62; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
% % plot3(points3(:,1),points3(:,2),points3(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
% % 
grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 6002; 
figure(fig_num); clf;
input_points = [points1; points2]; 
grid_size = 1;
grid_boundaries = [1 9 1 2 2.5 3.5]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));

XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 6002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(isempty(unmapped_grids))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%% 

fig_num = 789; 
figure(fig_num);clf;

[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(Velodyne_Lidar_Scan_ENU, fig_num); 

Max_x = ceil(Max_x); Max_y = ceil(Max_y); Max_z = ceil(Max_z);
Min_x = floor(Min_x); Min_y = floor(Min_y); Min_z = floor(Min_z);

% std(Velodyne_Lidar_Scan_ENU(:,3))

fig_num = 7002; 
figure(fig_num); clf;
input_points = Velodyne_Lidar_Scan_ENU; 
grid_size = 1;
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 
point_density = 10;
theta_threshold = pi/9; 
std_threshold = 0.02; 
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));

%%

XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

XYZ_matrix = unique(XYZ_matrix,'rows'); 

[~,XYZ_matrix_indices] = unique(XYZ_matrix(:,1:2),'rows'); 

XYZ_matrix = XYZ_matrix(XYZ_matrix_indices,:);  

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
% x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
% y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

x_range = unique(XYZ_matrix(:,1))'; 
y_range = unique(XYZ_matrix(:,2))';
% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% C = setdiff(XY_pairs, XY_pairs(idx==0,:), 'rows');

% Fill Z matrix with corresponding Z values
Z(idx~=0) = XYZ_matrix(idx(idx~=0), 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 7002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)==zeros(0,1))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)
