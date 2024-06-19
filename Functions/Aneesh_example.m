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
boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
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
boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
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
boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
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
boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary circle
angles = linspace(0,360,100)'*pi/180;
plot(radius*cos(angles),radius*sin(angles),'b-','LineWidth',3);

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
function boundary_points = fcn_INTERNAL_findBoundaryPoints(X,Y,Z,x_range,y_range, fig_num)
% Find the X_interval
x_interval = x_range(2)-x_range(1);
y_interval = y_range(2)-y_range(1);

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
    ylim([min(y_range) max(y_range)]);

    % Make axis slightly bigger than range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);


    % Plot the results
    plot(boundary_points_falling_y(:,1),boundary_points_falling_y(:,2),'r.','Markersize',40);
    plot(boundary_points_rising_y(:,1),boundary_points_rising_y(:,2),'g.','Markersize',40);
    plot(boundary_points_falling_x(:,1),boundary_points_falling_x(:,2),'c.','Markersize',40);
    plot(boundary_points_rising_x(:,1),boundary_points_rising_x(:,2),'m.','Markersize',40);
end


end % Ends fcn_INTERNAL_findBoundaryPoints
