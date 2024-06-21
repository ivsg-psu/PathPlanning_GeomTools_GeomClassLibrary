%% script_test_fcn_geometry_findBoundaryPoints
% - this is is a script written to test the function:
% fcn_geometry_findBoundaryPoints.m
%

% Revision history:
% 2024_06_21 - jpz5469@psu.edu
% -- original write of the code
close all;

%% Example for Aneesh

% Create some data
N_points = 10;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,15);

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));
% Make all Z data which has XY data above line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = Y>Y_line;
Z(flag_larger_than) = 1;
if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

% Calculate boundary points
fig_num = 1;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = x_range+2;
plot(x_range,y_line,'b-','LineWidth',3);

%% EXAMPLE 2
% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which has XY data below line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = find(Y<Y_line);
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

% Calculate boundary points
fig_num = 2;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = x_range+2;
plot(x_range,y_line,'b-','LineWidth',3);

%% EXAMPLE 3
% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which are near to origin equal to 1
distance_squared = X.^2 + Y.^2;
radius = 2;
radius_squared = radius.^2;

flag_larger_than = find(distance_squared<radius_squared);
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

% Calculate boundary points
fig_num = 3;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary circle
angles = linspace(0,360,100)'*pi/180;
plot(radius*cos(angles),radius*sin(angles),'b-','LineWidth',3);

%% EXAMPLE 4
% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,2,N_points);

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which are near to origin equal to 1
distance_squared = X.^2 + Y.^2;
radius = 2.5;
radius_squared = radius.^2;

flag_larger_than = find(distance_squared<radius_squared);
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

% Calculate boundary points
fig_num = 4;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary circle
angles = linspace(0,360,100)'*pi/180;
plot(radius*cos(angles),radius*sin(angles),'b-','LineWidth',3);
