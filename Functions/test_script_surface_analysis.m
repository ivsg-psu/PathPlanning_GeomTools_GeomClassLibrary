
% 2024_06_12 - Aneesh Batchu

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


pointDensity = 20;
[drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, pointDensity, (1));

flag_do_debug = 1; 

[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

% Calculate the number of points in each grid
total_N_points_in_each_grid = fcn_geometry_findRepeatedIndices(gridIndices,length(grid_AABBs(:,1)), -1);

% Save all the indices to a cell array. Each cell contains the indices of
% the inputPoints that belong to the grid. 
gridIndices_cell_array = fcn_geometry_createGridPointIndicesCellArray(gridIndices,length(grid_AABBs(:,1)),-1);

% Unmapped/not_fitted grids. These grids do not contain enough number of
% points to fit a plane and classify it as drivable or non-drivable
original_unmapped_grid_numbers = find(total_N_points_in_each_grid(:,1) < pointDensity); 

% Mapped grids. Later these grids are classified into drivable and non
% drivable
original_mapped_grid_numbers = find(total_N_points_in_each_grid(:,1) >= pointDensity); 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_grids = gridIndices_cell_array(original_mapped_grid_numbers); 

% Indices of points in mapped grids
indices_original_mapped_grids = vertcat(original_mapped_grids{:}); 

% Input points in the mapped grids
points_in_original_mapped_grids = inputPoints(indices_original_mapped_grids,:); 


if flag_do_debug

    plot_3D = 0; 

    fig_num = 100011; 
    figure(fig_num); clf;

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    % The indices of the mapped grids are extracted and concatenated
    original_unmapped_grids = gridIndices_cell_array(original_unmapped_grid_numbers);

    % Indices of points in mapped grids
    indices_original_unmapped_grids = vertcat(original_unmapped_grids{:});

    % Input points in the mapped grids
    points_in_original_unmapped_grids = inputPoints(indices_original_unmapped_grids,:);
 
    if plot_3D
        view(3)
        % Plot the unmapped points red
        plot3(points_in_original_unmapped_grids(:,1),points_in_original_unmapped_grids(:,2),points_in_original_unmapped_grids(:,3),'.','MarkerSize',20,'Color',[0.6350 0.0780 0.1840]);
        % Plot the mapped points green
        plot3(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),points_in_original_mapped_grids(:,3),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);
    else
        % Plot the unmapped points red
        plot(points_in_original_unmapped_grids(:,1),points_in_original_unmapped_grids(:,2),'.','MarkerSize',20,'Color',[0.6350 0.0780 0.1840]);
        % Plot the mapped points green
        plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);
    end


    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
end

% [Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(points_in_original_mapped_grids, (-1));
% 
% gridBoundaries = [Min_x,Max_x,Min_y,Max_y,Min_z,Max_z];
% 
% gridSize_mapped = round(min([Max_x - Min_x, Max_y - Min_y, Max_z - Min_z])/2,2); 
% 
% [gridIndices_mapped,grid_AABBs,~] = fcn_geometry_separatePointsIntoGrids(points_in_original_mapped_grids, gridSize_mapped, gridBoundaries, (fig_num+1));
% 
% % Calculate the number of points in each grid
% total_points_in_each_grid = fcn_geometry_findRepeatedIndices(gridIndices_mapped,length(grid_AABBs(:,1)), -1);
% 
% % Save all the indices to a cell array. Each cell contains the indices of
% % the inputPoints that belong to the grid. 
% gridIndices_cell_array_mapped = fcn_geometry_createGridPointIndicesCellArray(gridIndices_mapped,length(grid_AABBs(:,1)),-1);

total_mapped_grids = length(original_mapped_grids); 

parameters_of_fitted_plane = zeros(total_mapped_grids,3); 

unit_normal_vectors = zeros(total_mapped_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 


for ith_mapped_grid = 1:total_mapped_grids
    [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_mapped_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_mapped_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(original_mapped_grids{ith_mapped_grid},:),-1);
        parameters_of_fitted_plane(ith_mapped_grid,:) = parameters';
end

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < 0.1; 

% Grids that satisy both the conditions (STEP 1 and STEP 2). The grids that
% are within the vertical and std threshold
mapped_grids_within_vertical_and_std_thresholds = (mapped_grids_within_vertical_threshold == 1) & (mapped_grids_within_std_threshold == 1);

% Find the drivable grids (original)
drivable_grids = original_mapped_grid_numbers(mapped_grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids (original)
non_drivable_grids = original_mapped_grid_numbers(mapped_grids_within_vertical_and_std_thresholds == 0);

% Final drivable grid numbers of the mapped grids
drivable_grid_numbers_in_mapped_grids = find(ismember(original_mapped_grid_numbers, drivable_grids));

% Final non drivable grid numbers of the mapped grids
non_drivable_grid_numbers_in_mapped_grids = find(ismember(original_mapped_grid_numbers, non_drivable_grids));


% % Plot all the points
% plot(inputPoints(:,1),inputPoints(:,2),'.','MarkerSize',20,'Color',[0 0 0]);

if flag_do_debug
    fig_num = 100012;
    figure(fig_num); clf;

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    view(3)

    % Plot all the points
    plot3(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),points_in_original_mapped_grids(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

    for ith_domain = 1:length(drivable_grids)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
        current_color = [0.4660 0.6740 0.1880];

        % Plot current AABB
        current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

        % Nudge the current AABB inward
        current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

        % Calculate the gridlines
        gridlines = [...
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan];

        % Plot the result
        plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


        % Get all points in this domain and plot them
        rows_in_domain = gridIndices==drivable_grids(ith_domain);
        points_in_domain = inputPoints(rows_in_domain,:);
        plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);

    % Plot the unit normal vector
    quiver3(gridCenters(drivable_grids(ith_domain),1),gridCenters(drivable_grids(ith_domain),2),gridCenters(drivable_grids(ith_domain),3), unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),1),unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),2),unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),3),0,'g','Linewidth',3);
    
    end

    for ith_domain = 1:length(non_drivable_grids)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
        current_color = [0.6350 0.0780 0.1840];

        % Plot current AABB
        current_AABB = grid_AABBs(non_drivable_grids(ith_domain),:);

        % Nudge the current AABB inward
        current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

        % Calculate the gridlines
        gridlines = [...
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
            nan nan nan;
            current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
            current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
            nan nan nan];

        % Plot the result
        plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


        % Get all points in this domain and plot them
        rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
        points_in_domain = inputPoints(rows_in_domain,:);
        plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);

    % Plot the unit vector
    quiver3(gridCenters(non_drivable_grids(ith_domain),1),gridCenters(non_drivable_grids(ith_domain),2),gridCenters(non_drivable_grids(ith_domain),3), unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),1),unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),2),unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),3),0,'g','Linewidth',3);
    end

    % % Plot the unit vector
    % quiver3(base_point(1,1),base_point(1,2),base_point(1,3), unit_vector(1,1),unit_vector(1,2),unit_vector(1,3),0,'g','Linewidth',3);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
end

if flag_do_debug
    fig_num = 100012;
    figure(fig_num); clf; 

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0, 0, 0]);

% Plot the input points by domain with different colors for each
% domain
for ith_domain = 1:length(drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.4660 0.6740 0.1880];
    % Plot current AABB
    current_AABB = grid_AABBs(length(drivable_grids),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Plot the result
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',20,'Color',current_color);
end

% Plot the input points by domain with different colors for each
% domain
for ith_domain = 1:length(non_drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.6350 0.0780 0.1840];
    % Plot current AABB
    current_AABB = grid_AABBs(length(non_drivable_grids),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Plot the result
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',20,'Color',current_color);
end
% Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
end

% [drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, (fig_num)); 


%% Different data - drivable and non_drivable case

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
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points3(:,1),points3(:,2),points3(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 333; 
figure(fig_num); clf;
inputPoints = [points1; points2; points3]; 
gridSize = 1;
gridBoundaries = [1 4 1 2 2.5 5.5]; 

[drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, (123)); 


[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

% Counter to the grids with more than 20 points
grids_with_sufficient_points= zeros(total_N_grids,1); 

% Points per grid (point density), think about delta Y
point_density = 20; 


for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));

    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= point_density
        grids_with_sufficient_points(ith_grid,1) = 1; 
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        
    end
end

% Unmapped/not fitted grids 
unmapped_grids = find(total_N_points_in_each_grid(:,1) < point_density);

% Mapped/Fitted grids 
mapped_grids = find(total_N_points_in_each_grid(:,1) >= point_density);

% Find the grids with a fitted plane. Simply, a grid with more than 20
% points
% grids_with_fitted_plane = find(grids_with_sufficient_points == 1);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors(mapped_grids,:).*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
grids_within_std_threshold = standard_deviation_in_plane_orthogonals(mapped_grids,1) < 0.1; 

% Grids that satisy both the conditions (STEP 1 and STEP 2). The grids that
% are within the vertical and std threshold
grids_within_vertical_and_std_thresholds = (grids_within_vertical_threshold == 1) & (grids_within_std_threshold == 1);

% Find the drivable grids
drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids
non_drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds == 0);

% Detect edges using a method like Sobel
edges = edge(grids_within_vertical_and_std_thresholds, 'Sobel');

% Visualize the edges on top of the original grid map
figure(1);
imshow(edges);
title('Detected Boundaries between Drivable and Non-Drivable Surfaces');


fig_num = 90909; 
figure(fig_num);clf;
hold on
view(3)
plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
for ith_domain = 1:length(drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.4660 0.6740 0.1880]; 

    % Plot current AABB
    current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(non_drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.6350 0.0780 0.1840]; 

    % Plot current AABB
    current_AABB = grid_AABBs(non_drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(unmapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.5 0.5 0.5];

    % Plot current AABB
    current_AABB = grid_AABBs(unmapped_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==unmapped_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end


 %%



N_points = 300;
Ext_Square_Size=10;
Int_Square_Size=3;
fig_num = 31;
diag_flag=1;
noise= 0.15;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);

% Surface Analysis

% fig_num = 333; 
% figure(fig_num); clf;
inputPoints = points; 
gridSize = 1;
gridBoundaries = [-5 5 -5 5 -1 1]; 

[drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, (123)); 


%% Create data

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

true_parameters = [ 1 0 1]';
points = [2 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 

fig_num = 1; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 333; 
figure(fig_num); clf;
inputPoints = [points1; points2]; 
gridSize = 1;
gridBoundaries = [1 3 1 2 3 4]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));


    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= 20
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        % [parameters, standard_deviation_in_z, z_fit, unit_normal_vector(ith_grid), base_point, standard_deviation_in_plane_orthogonals(ith_grid)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);

    end
end

% Find the grids with a fitted plane
grids_with_fitted_plane = find(standard_deviation_in_plane_orthogonals ~= 0);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = find(angle_btw_unit_normals_and_vertical < pi/9);

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within vertical and standard deviation limit
drivable_grids = find(standard_deviation_in_plane_orthogonals(grids_within_vertical_threshold,:) < 0.1); 
            


%% 900 points: all drivable surfaces

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

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

% Counter to the grids with more than 20 points
grids_with_sufficient_points= zeros(total_N_grids,1); 

% Points per grid (point density), think about delta Y
point_density = 20; 


for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));

    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= point_density
        grids_with_sufficient_points(ith_grid,1) = 1; 
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        
    end
end

% Unmapped/not fitted grids 
unmapped_grids = find(total_N_points_in_each_grid(:,1) < point_density);

% Mapped/Fitted grids 
mapped_grids = find(total_N_points_in_each_grid(:,1) >= point_density);

% Find the grids with a fitted plane. Simply, a grid with more than 20
% points
% grids_with_fitted_plane = find(grids_with_sufficient_points == 1);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors(mapped_grids,:).*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
grids_within_std_threshold = standard_deviation_in_plane_orthogonals(mapped_grids,1) < 0.1; 

% Grids that satisy both the conditions (STEP 1 and STEP 2). The grids that
% are within the vertical and std threshold
grids_within_vertical_and_std_thresholds = (grids_within_vertical_threshold == 1) & (grids_within_std_threshold == 1);

% Find the drivable grids
drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids
non_drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds == 0);

% Detect edges using a method like Sobel
edges = edge(grids_within_vertical_and_std_thresholds, 'Sobel');


% set unmapped grids to -1 (might use this to plot or in any other
% operation)
% grids_with_sufficient_points(unmapped_grids,1) = -1;  


% Final Plotting

fig_num = 90909; 
figure(fig_num);clf;
hold on
view(3)
plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
for ith_domain = 1:length(drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.4660 0.6740 0.1880]; 

    % Plot current AABB
    current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(non_drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.6350 0.0780 0.1840]; 

    % Plot current AABB
    current_AABB = grid_AABBs(non_drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(unmapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.5 0.5 0.5];

    % Plot current AABB
    current_AABB = grid_AABBs(unmapped_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==unmapped_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end



%% Different data - 

% Create data

rng(123)

N_points = 100; 

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

N_points = 100; 

true_parameters = [ 0 0 3]';
points = [5 1].*ones(N_points,1)+ [4,1].*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

% Add noise
true_sigma = 0.15;
z = true_z + true_sigma*randn(N_points,1);

points2 = [x, y, z]; 


fig_num = 1; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 1]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 333; 
figure(fig_num); clf;
inputPoints = [points1; points2; points3]; 
gridSize = 1;
gridBoundaries = [1 4 1 2 2.5 5.5]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

% Counter to the grids with more than 20 points
grids_with_sufficient_points= zeros(total_N_grids,1); 

% Points per grid (point density), think about delta Y
point_density = 20; 


for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));

    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= point_density
        grids_with_sufficient_points(ith_grid,1) = 1; 
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        
    end
end

% Unmapped/not fitted grids 
unmapped_grids = find(total_N_points_in_each_grid(:,1) < point_density);

% Mapped/Fitted grids 
mapped_grids = find(total_N_points_in_each_grid(:,1) >= point_density);

% Find the grids with a fitted plane. Simply, a grid with more than 20
% points
% grids_with_fitted_plane = find(grids_with_sufficient_points == 1);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors(mapped_grids,:).*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
grids_within_std_threshold = standard_deviation_in_plane_orthogonals(mapped_grids,1) < 0.1; 

% Grids that satisy both the conditions (STEP 1 and STEP 2). The grids that
% are within the vertical and std threshold
grids_within_vertical_and_std_thresholds = (grids_within_vertical_threshold == 1) & (grids_within_std_threshold == 1);

% Find the drivable grids
drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids
non_drivable_grids = mapped_grids(grids_within_vertical_and_std_thresholds == 0);


fig_num = 90909; 
figure(fig_num);clf;
hold on
view(3)
plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
for ith_domain = 1:length(drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.4660 0.6740 0.1880]; 

    % Plot current AABB
    current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(non_drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.6350 0.0780 0.1840]; 

    % Plot current AABB
    current_AABB = grid_AABBs(non_drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end

for ith_domain = 1:length(unmapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.5 0.5 0.5];

    % Plot current AABB
    current_AABB = grid_AABBs(unmapped_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==unmapped_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end


%% Create data

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

true_parameters = [ 1 0 1]';
points = [2 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 

fig_num = 1; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 333; 
figure(fig_num); clf;
inputPoints = [points1; points2]; 
gridSize = 1;
gridBoundaries = [1 3 1 2 3 4]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));


    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= 20
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        % [parameters, standard_deviation_in_z, z_fit, unit_normal_vector(ith_grid), base_point, standard_deviation_in_plane_orthogonals(ith_grid)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);

    end
end

% Find the grids with a fitted plane
grids_with_fitted_plane = find(standard_deviation_in_plane_orthogonals ~= 0);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = find(angle_btw_unit_normals_and_vertical < pi/9);

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within vertical and standard deviation limit
drivable_grids = find(standard_deviation_in_plane_orthogonals(grids_within_vertical_threshold,:) < 0.1); 
            


%% 900 points: all drivable surfaces

rng(123)

N_points = 1500; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ 9*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points1 = [x, y, z]; 

fig_num = 1011;
figure(fig_num)

view(3)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')


fig_num = 10333; 
figure(fig_num); clf;
inputPoints = points1;
gridSize = 1;
gridBoundaries = [1 10 1 10 2.5 3.5]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);

parameters_of_fitted_plane = zeros(total_N_grids,3); 

unit_normal_vectors = zeros(total_N_grids,3); 

standard_deviation_in_plane_orthogonals = zeros(total_N_grids,1); 

% Counter to the grids with more than 20 points
grids_with_sufficient_points= zeros(total_N_grids,1); 

% Points per grid (point density)
point_density = 20; 

for ith_grid = 1:total_N_grids

    % Find the no.of points in each grid
    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid));

    % Fit a lane to a grid that contains more than 20 points
    if total_N_points_in_each_grid(ith_grid,1) >= point_density
        grids_with_sufficient_points(ith_grid,1) = 1; 
        [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
        parameters_of_fitted_plane(ith_grid,:) = parameters';
        
    end
end

% Find the grids with a fitted plane. Simply, a grid with more than 20
% points
grids_with_fitted_plane = find(grids_with_sufficient_points == 1);  

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
grids_within_std_threshold = standard_deviation_in_plane_orthogonals < 0.1; 

% Find the drivable grid number
drivable_grids = find(grids_within_vertical_threshold==grids_within_std_threshold); 

% Unmapped/not fitted grids 
unmapped_grids = find(total_N_points_in_each_grid(:,1) < point_density);

% set unmapped grids to -1 (might use this to plot or in any other
% operation)
grids_with_sufficient_points(unmapped_grids,1) = -1;  


% Final Plotting

fig_num = 90909; 
figure(fig_num);clf;
hold on
view(3)
plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
for ith_domain = 1:length(drivable_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.4660 0.6740 0.1880]; 

    % Plot current AABB
    current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

    % Nudge the current AABB inward
    current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
        nan nan nan;
        current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
        current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
        nan nan nan];

    % Plot the result
    plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==drivable_grids(ith_domain);
    points_in_domain = inputPoints(rows_in_domain,:);
    plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
end



%% Create data

N_points = 100; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points1 = [x, y, z]; 

true_parameters = [ 1 0 1]';
points = [2 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 

fig_num = 1; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')


%% Create data

% Define true parameters
% true_parameters = [0; 0; 1]; % Change the z coefficient to 1

% Generate points
points = 1*ones(100,2) + rand(100,2);
x = points(:,1);
y = points(:,2);

% Generate z values within the range [1, 2]
true_z = linspace(1,2,100)'; % Generate random z values in the range [1, 2]

% Create the points matrix
points1 = [x, y, true_z];



% Plot the points
figure(1)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

%% Create data

true_parameters = [ 0 0 2]';
points = 1*ones(100,2) + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points1 = [x, y, z]; 


true_parameters = [ 0 0 3]';
points = [ones(100,1), 2*ones(100,1)] + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 


plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

grid on

xlabel('x')
ylabel('y')
zlabel('z')

%% BASIC example 1: XY with lots of data
fig_num = 200;
figure(fig_num);
clf;

inputPoints = [10*rand(300,2), ones(300,1)];
gridSize = 1;
gridBoundaries = [1 5 1 5 0 2]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


%% Seperate data into grids

% Create data

true_parameters = [ 0 0 0]';
points = 1*ones(100,2) + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points1 = [x, y, z]; 


true_parameters = [ 0 0 1]';
points = [ones(100,1), 2*ones(100,1)] + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 


fig_num = 333; 
inputPoints = [points1(:,1:2); points2(:,1:2)]; 
gridSize = 0.5;
gridBoundaries = [1 2 1 3]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

%% Seperate data into grids

% Create data

true_parameters = [ 0 0 0]';
points = 1*ones(100,2) + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points1 = [x, y, z]; 


true_parameters = [ 0 0 1]';
points = [ones(100,1), 2*ones(100,1)] + rand(100,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points2 = [x, y, z]; 


fig_num = 333; 
inputPoints = [points1; points2]; 
gridSize = 1;
gridBoundaries = [1 2 1 3 0 2]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

total_N_grids = length(grid_AABBs(:,1));

total_N_points_in_each_grid = zeros(total_N_grids,1);
for ith_grid = 1:total_N_grids

    points_in_each_grid = find(gridIndices==ith_grid);
    total_N_points_in_each_grid(ith_grid,1) = length(find(gridIndices==ith_grid)); 
    
    if length(find(gridIndices==ith_grid)) >= 20
        [parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_fitPlaneLinearRegression(inputPoints(points_in_each_grid,:),fig_num);
    end
end




%% Test 2: add noise

fig_num = 1;
figure(fig_num);
clf;
rng(1823);

N_points = 100; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

true_sigma = 0.1;
z = true_z + true_sigma*randn(N_points,1);

% z = true_z; 

points1 = [x, y, z]; 

fig_num = 22; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')


fig_num = 1;
figure(fig_num);
clf;
[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_fitPlaneLinearRegression([x y z],fig_num);

% assert(isequal(round(true_parameters,1),round(fitted_parameters,1)));
% assert(isequal(round(true_sigma,1),round(standard_deviation_in_z,1)));


