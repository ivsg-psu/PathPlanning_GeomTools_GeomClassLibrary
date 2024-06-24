
close all; 

% Revision History
% Aneesh Batchu - 2024_06_21
% -- wrote the code originally

%% Basic Test: grid_size = 1. Less data, all the surfaces are mapped: All the surfaces have enough data to fit a plane

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

fig_num = 1; clf;
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
x_range = min(gridCenters_mapped_grids(:,1)):grid_size:max(gridCenters_mapped_grids(:,1)); % min_x:gridSize:max_x 
y_range = min(gridCenters_mapped_grids(:,2)):grid_size:max(gridCenters_mapped_grids(:,2)); % min_y:gridSize:max_y 

% % Given x_range and y_range
% x_range = min(gridCenters_mapped_grids(:,1)):grid_size:max(gridCenters_mapped_grids(:,1)); % min_x:gridSize:max_x 
% y_range = min(gridCenters_mapped_grids(:,2)):grid_size:max(gridCenters_mapped_grids(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = zeros(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% Display the result
disp('Z = ');
disp(Z);

x_limits = [min(x_range) max(x_range)];  
y_limits = [0 2]; 
% Calculate boundary points
fig_num = 1002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

% boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range,grid_size,fig_num);
% plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

%% 

rng (123)

N_points = 300;
Ext_Square_Size=3;
Int_Square_Size=1;
fig_num = 32;
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
x_range = min(gridCenters_mapped_grids(:,1)):grid_size:max(gridCenters_mapped_grids(:,1)); % min_x:gridSize:max_x 
y_range = min(gridCenters_mapped_grids(:,2)):grid_size:max(gridCenters_mapped_grids(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = zeros(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% Display the result
disp('Z = ');
disp(Z);

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 
x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 2002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

%% 

rng (123)

N_points = 2500;
Ext_Square_Size=4;
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
grid_size = 0.4;
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
x_range = [-0.6000   -0.2000    0.2000    0.6000]; min(gridCenters_mapped_grids(:,1)):grid_size:max(gridCenters_mapped_grids(:,1)); % min_x:gridSize:max_x 
y_range = [-0.6000   -0.2000    0.2000    0.6000];%min(gridCenters_mapped_grids(:,2)):grid_size:max(gridCenters_mapped_grids(:,2)); % min_y:gridSize:max_y 

% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = zeros(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% Fill Z matrix with corresponding Z values
Z(:) = XYZ_matrix(idx, 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% Display the result
disp('Z = ');
disp(Z);

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
fig_num = 3002;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

%% Basic Test: grid_size = 0.5. Similar to the previous case. The drivable and non-drivable 

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

fig_num = 1; clf;
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

fig_num = 2002; 
figure(fig_num); clf;
input_points = [points1; points2; points3]; 
grid_size = 0.2;
grid_boundaries = [1 4 1 2 2.5 5.5]; 
point_density = 20;
theta_threshold = pi/9; 
std_threshold = 0.1; 
flag_plot_in_3D = 0; 

[drivable_grids, non_drivable_grids, unmapped_grids, gridCenters_mapped_grids, drivable_grid_numbers_in_mapped_grids, non_drivable_grid_numbers_in_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, (flag_plot_in_3D), (fig_num));


%% 

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

fig_num = 336; 
figure(fig_num); clf;
inputPoints = points; 
gridSize = 1;
gridBoundaries = [-3 3 -3 3 -1 1]; 

[drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, (fig_num)); 


%% Different data - drivable and non_drivable case

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

fig_num = 3; clf;
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

fig_num = 335; 
figure(fig_num); clf;
inputPoints = [points1; points2]; 
gridSize = 1;
gridBoundaries = [1 9 1 2 2.5 3.5]; 

[drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, (fig_num)); 




%% Functions start here
function [boundary_points_falling, boundary_points_rising] = fcn_INTERNAL_findBoundaryPointsX(X, Y, Z, y_interval)

flag_do_debug = 0;

Z_greater_than = find(Z>0.5);

% Pad the first and last points
N_points = numel(X); % Total number of elements
Z_greater_than_padded = [0; Z_greater_than; N_points+1];

%% FInd changes
changes_in_sequence = diff(Z_greater_than_padded);
indicies_falling_edge = find(changes_in_sequence>1.5); %  Anything greater than 1 is a change. Have to add 1 because we padded indicies above
indicies_rising_edge = find(changes_in_sequence>1.5)+1; % Anything greater than 1 is a change. Have to add 1 because we padded indicies above

indicies_with_rising_x_edge  = Z_greater_than_padded(indicies_rising_edge);
indicies_with_falling_x_edge = Z_greater_than_padded(indicies_falling_edge);

% Clean up any indicies outside of range
indicies_with_rising_x_edge(indicies_with_rising_x_edge<1) = [];
indicies_with_rising_x_edge(indicies_with_rising_x_edge>N_points) = [];
indicies_with_falling_x_edge(indicies_with_falling_x_edge<1) = [];
indicies_with_falling_x_edge(indicies_with_falling_x_edge>N_points) = [];

% Clean up any indicies on borders
border_dimension = length(X(:,1));
border_indices_rising = find(mod(indicies_with_rising_x_edge,border_dimension)==1);
indicies_with_rising_x_edge(border_indices_rising) = [];
border_indices_falling = find(mod(indicies_with_falling_x_edge,border_dimension)==0);
indicies_with_falling_x_edge(border_indices_falling) = [];




boundary_points_rising  = [X(indicies_with_rising_x_edge),Y(indicies_with_rising_x_edge)-y_interval/2];
boundary_points_falling = [X(indicies_with_falling_x_edge),Y(indicies_with_falling_x_edge)+y_interval/2];

% Plot the data in 2D?
if 1==flag_do_debug
    figure(1111);
    clf;
    hold on;

    % Plot the inputs
    plot(X(1:N_points),Y(1:N_points),'.','Color',[0.5 0.5 0.5],'Markersize',20);
    plot(X(Z_greater_than),Y(Z_greater_than),'k.','Markersize',50);

    % Number the results (for clarity)
    for ith_label = 1:length(Z_greater_than)
        label_number = Z_greater_than(ith_label);
        current_text = sprintf('%.0d',label_number);
        text(X(label_number),Y(label_number),current_text,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');
    end
    xlabel('X [m]');
    ylabel('Y [m]');

    % Start by forcing tight axes
    xlim([min(X,[],'all') max(X,[],'all')]);
    ylim([min(Y,[],'all') max(Y,[],'all')]);

    % Make axis slightly bigger than range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    plot(boundary_points_falling(:,1),boundary_points_falling(:,2),'r.','Markersize',30);
    plot(boundary_points_rising(:,1),boundary_points_rising(:,2),'g.','Markersize',30);
end

end

%% fcn_INTERNAL_findBoundaryPoints
function boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range,grid_size, fig_num)
% Find the X_interval
% x_interval = x_range(2)-x_range(1);
% y_interval = y_range(2)-y_range(1);

x_interval = grid_size;
y_interval = grid_size; 

[boundary_points_falling_y, boundary_points_rising_y] = fcn_INTERNAL_findBoundaryPointsX(X, Y, Z,  y_interval);
[boundary_points_falling_x_transpose, boundary_points_rising_x_transpose] = fcn_INTERNAL_findBoundaryPointsX(Y', X', Z', x_interval);
boundary_points_falling_x = fliplr(boundary_points_falling_x_transpose);
boundary_points_rising_x = fliplr(boundary_points_rising_x_transpose);

boundary_points = [boundary_points_falling_y; boundary_points_rising_y; boundary_points_falling_x; boundary_points_rising_x];

if 1 == 1
    % Plot the data in 2D
    figure(fig_num);
    clf;
    hold on;
    axis equal;

    % Plot the results
    flag_larger_than = Z>0.5;
    plot(X(flag_larger_than),Y(flag_larger_than),'k.','Markersize',50);

    xlabel('X [m]');
    ylabel('Y [m]');

    xlim([min(x_range) max(x_range)]);
    % ylim([min(y_range) max(y_range)]);
    ylim([min(y_range) max(y_range)])
    % ylim([0 2]);

    % Make axis slightly bigger than range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    if isempty(boundary_points_falling_y)
        boundary_points_falling_y = zeros(0,2);
    end
    if isempty(boundary_points_rising_y)
        boundary_points_rising_y = zeros(0,2);
    end
    if isempty(boundary_points_falling_x)
        boundary_points_falling_x = zeros(0,2);
    end
    if isempty(boundary_points_rising_x)
        boundary_points_rising_x = zeros(0,2);
    end

    % Plot the results
    plot(boundary_points_falling_y(:,1),boundary_points_falling_y(:,2),'r.','Markersize',40);
    plot(boundary_points_rising_y(:,1),boundary_points_rising_y(:,2),'g.','Markersize',40);
    plot(boundary_points_falling_x(:,1),boundary_points_falling_x(:,2),'c.','Markersize',40);
    plot(boundary_points_rising_x(:,1),boundary_points_rising_x(:,2),'m.','Markersize',40);
end


end % Ends fcn_INTERNAL_findBoundaryPoints

