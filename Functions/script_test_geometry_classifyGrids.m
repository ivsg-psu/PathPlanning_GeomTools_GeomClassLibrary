

%% INPUTS

% figure number
fig_num_LLA = 3001;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% As of now. 
point_density = 60;

%% MAIN - grids with zero and non-zero

% INPUTS - input_points,grid_size,grid_boundaries, plotting
% OUTPUTS - gridIndices_cell_array, 
% total_N_points_in_each_grid,gridCenters, grids_with_zero_points, grids_greater_than_zero_points 
% fcn_geometry_findGridsWithPoints

% Revision History
% 2024_06_15 - Aneesh Batchu
% -- Wrote the code originally
% 2024_06_25 - Aneesh Batchu
% -- NaNs from gridIndices are removed before they are used in the surface
% analysis
% 2024_07_15 - Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis
% 2024_07_15 - Jiabao 
% -- Functionalized this code

% Divides the data into grids
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(input_points, grid_size, grid_boundaries, (-1));

% Remove NaNs from gridIndices if there any
gridIndices = gridIndices(~isnan(gridIndices));

% Save all the indices to a cell array. Each cell contains the indices of
% the inputPoints that belong to the grid. 
% The second input length(gridCenters(:,1)) is to assign the size of an array. Eaxh grid has a grid center.  
gridIndices_cell_array = fcn_geometry_createGridPointIndicesCellArray(gridIndices,length(gridCenters(:,1)),-1);

% Calculate the number of points in each grid
total_N_points_in_each_grid = fcn_geometry_findRepeatedIndices(gridIndices,length(grid_AABBs(:,1)), -1);

% Find grids with zero points (empty grid numbers)
grids_with_zero_points = find((total_N_points_in_each_grid(:,1) == 0)); 

% Grid Centers of the grids with zero point density 
gridCenters_zero_point_density = gridCenters(grids_with_zero_points,1:3); 

% Find grids with more than zero points
grids_greater_than_zero_points = find((total_N_points_in_each_grid(:,1) > 0)); 

% Grid Centers of the grids with zero point density (Unmapped grid centers)
gridCenters_greater_than_zero_point_density = gridCenters(grids_greater_than_zero_points,1:3); 

% Find the point density based on the Histogram analysis
% [point_density,actual_driving_surface_grids_hist,total_grids_hist] = fcn_geometry_findPointDensity(total_points_in_each_grid_in_actual_driving_surface,total_points_in_each_grid,fig_num)


% plotting
 
fig_num = 11; 
figure(fig_num);clf

marker_size = 10;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with zero point density';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_zero_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids greater than zero point density';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_greater_than_zero_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);


%% MAIN - classify grids with more than zero points into mapped and unmapped


% INPUTS -
% point_density,total_N_points_in_each_grid,grids_greater_than_zero_points,
% grid_centers, plotting options
% OUTPUTS - original_grids_with_required_point_density
% fcn_geometry_classifyGridsintoMappedUnmapped

% Revision History
% 2024_06_15 - Aneesh Batchu
% -- Wrote the code originally
% 2024_07_15 - Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis

% Find grids with low point density but not zero point density (Unmapped grid centers)
% Original: The grid numbers belong to the initial grid indices cell array
% The numbering does not start from "1"
original_grids_with_low_point_density = grids_greater_than_zero_points((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Find mapped grids 
% Mapped grids. Later these grids are classified into drivable and non
% drivable
% Original: The grid numbers belong to the initial grid indices cell array
% The numbering does not start from "1"
original_grids_with_required_point_density = grids_greater_than_zero_points(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% Grid Centers of the grids with low point density (Unmapped grid centers)
gridCenters_low_point_density = gridCenters(original_grids_with_low_point_density,1:3);  

% Grid Centers of the grids with required point density (Mapped grid centers)
gridCenters_required_point_density = gridCenters(original_grids_with_required_point_density,1:3);  

% Current grid numbers of the grids with low point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_low_point_density = find((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Current grid numbers of the grids with required point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_required_point_density = find(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% plotting
fig_num = 11; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Unmapped grids';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_low_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Mapped grids';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_required_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% ENU - plotting mapped and unmapped
figure(1234)
clf;
plot(gridCenters_low_point_density(:,1), gridCenters_low_point_density(:,2), '.','MarkerSize',40,'Color',[0.8 0.8 0.8]);
hold on
plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_grids_with_low_point_density(:,1))
    current_text = sprintf('%.0d',current_grids_with_low_point_density(ith_text));
    % Place the text on the grid center
    text(gridCenters_low_point_density(ith_text,1), gridCenters_low_point_density(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 8);
end

for ith_text = 1:length(current_grids_with_required_point_density(:,1))
    current_text = sprintf('%.0d',current_grids_with_required_point_density(ith_text));
     % Place the text on the grid center
    text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8);
end

%% MAIN - classify mapped grids into drivable and non-drivable

% INPUT
std_threshold = 0.1; 
theta_threshold = 30*pi/180;
z_diff_threshold = 0.1;

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_grids_with_required_point_density); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

% z_height of all the points 
mean_z_of_mapped_grids = zeros(total_mapped_grids,1); 

% Loop through all the mapped grids, recording standard deviation, unit
% vectors 
% if 0==flag_max_speed
%     h_waitbar = waitbar(0,'Performing surface analysis...');
% end

for ith_mapped_grid = 1:total_mapped_grids
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% STEP 1: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit)
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
% mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < std_threshold;
mapped_grids_within_std_threshold = standard_deviation_in_z < std_threshold;

% STEP 2: Find the mean of z height difference and classify the grid as
% drivable, if within the z_height_diff threshold.

% sort z values to find the first few max values
[sorted_z_diff_mapped_grids,sorted_indices_z_diff_mapped_grids] = sort(z_diff_mapped_grids,'descend');

% aa = sorted_indices_z_diff_mapped_grids(sorted_z_diff_mapped_grids>z_diff_threshold);
aa = (sorted_z_diff_mapped_grids>z_diff_threshold);

% STEP 3
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product);

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < theta_threshold;

% mapped_grids_within_vertical_and_std_thresholds = (mapped_grids_within_vertical_threshold == 1) & (mapped_grids_within_std_threshold == 1);
mapped_grids_within_vertical_and_std_thresholds = (aa == 0);

% Find the drivable grids (original)
drivable_grids = original_grids_with_required_point_density(mapped_grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids (original)
non_drivable_grids = original_grids_with_required_point_density(mapped_grids_within_vertical_and_std_thresholds == 0);

% Final drivable grid numbers of the mapped grids
drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, drivable_grids));

% Final non drivable grid numbers of the mapped grids
non_drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, non_drivable_grids));

% Grid centers of drivable grids 
gridCenters_drivable_grids = [gridCenters(drivable_grids,1), gridCenters(drivable_grids,2), gridCenters(drivable_grids,3), ones(length(drivable_grids),1)]; 

% Grid centers of nondrivable grids
gridCenters_non_drivable_grids = [gridCenters(non_drivable_grids,1), gridCenters(non_drivable_grids,2), gridCenters(non_drivable_grids,3), zeros(length(non_drivable_grids),1)]; 

% Concatenate the grid centers of drivable and non-drivable grids (2D)
gridCenters_mapped_grids = [gridCenters_drivable_grids; gridCenters_non_drivable_grids];

% plotting
fig_num = 11; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0, 1, 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_drivable_grids(:,1:3),marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [1, 0, 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
[~] = fcn_geometry_plotGridCenters(gridCenters_non_drivable_grids(:,1:3),marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% ENU mapped grids

figure(12345)
clf;
hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of the mapped grids')
plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1)); 
for ith_text = 1:length(gridCenters_required_point_density(:,1))
    current_text = sprintf('%.0d',mapped_grid_numbers(ith_text));
     % Place the text on the grid center
    text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8);
end

%% In-line scan z_height_difference

% figure number
fig_num_LLA = 3001;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1409});
% LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{:});
Trace_coordinates = LiDAR_outer_edge;

% Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_data = gps_object.ENU2WGSLLA(Trace_coordinates);

% Plot the ENU results
figure(fig_num_LLA);clf;

xlabel('Latitude')
ylabel('Longitude')

geoplot(LLA_data(:,1),LLA_data(:,2),'mo','MarkerSize',10);
hold on
geoplot(LLA_data(:,1),LLA_data(:,2),'k.','MarkerSize',10);

title('LLA Trace geometry')

geobasemap satellite
geotickformat -dd

%% Slide 4: LiDAR data is plotted ion ENU coordinates and z-coordinate is plotted in colors

fig_num_slide_4 = 4000; 

% LiDAR data
x = LiDAR_outer_edge(:,1);
y = LiDAR_outer_edge(:,2);
z = LiDAR_outer_edge(:,3);

% Height of the points
z_diff = abs(diff(z)); 

% sort z values to find the first few max values
[sorted_z_diff,sorted_indices_z_diff] = sort(z_diff,'descend'); 

% Normalize colors to the range [0, 1]
normalized_colors = (z_diff - min(z_diff)) / (max(z_diff) - min(z_diff));


aa = sorted_indices_z_diff(sorted_z_diff>0.06 & sorted_z_diff<0.07);


figure(fig_num_slide_4);clf;
% Create the scatter plot
scatter3(x(2:end), y(2:end), z(2:end), 20, normalized_colors, 'filled');
hold on
plot3(LiDAR_outer_edge(aa,1),LiDAR_outer_edge(aa,2),LiDAR_outer_edge(aa,3),'r.', 'MarkerSize',40)
colormap('jet'); 
colorbar; 
clim([0 1]); % Set color axis limits to match the normalized range
view(2); 
title('Z-coordinate is represented in colors');
xlabel('X[m]');
ylabel('Y[m]');


% figure number
fig_num_LLA = 3001;

% Plot the ENU results
figure(fig_num_LLA);


% geoplot(LLA_data(aa,1),LLA_data(aa,2),'c.','MarkerSize',20);
% hold on
% geoplot(LLA_data(aa,1),LLA_data(aa,2),'b.','MarkerSize',10);

geoplot(LLA_data(aa(end),1),LLA_data(aa(end),2),'c.','MarkerSize',20);
hold on
geoplot(LLA_data(aa(end),1),LLA_data(aa(end),2),'b.','MarkerSize',10);


%% 

LiDAR_scans = 1400:1404; 

% Total number of mapped grids
total_LiDAR_scans = length(LiDAR_scans);

for ith_LiDAR_scan = 1:total_LiDAR_scans
    LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{LiDAR_scans(ith_LiDAR_scan)});

    % LiDAR data
    x = LiDAR_outer_edge(:,1);
    y = LiDAR_outer_edge(:,2);
    z = LiDAR_outer_edge(:,3);

    % Height of the points
    z_diff = abs(diff(z));

    % sort z values to find the first few max values
    [sorted_z_diff,sorted_indices_z_diff] = sort(z_diff,'descend');

    % Normalize colors to the range [0, 1]
    normalized_colors = (z_diff - min(z_diff)) / (max(z_diff) - min(z_diff));


    aa = sorted_indices_z_diff(sorted_z_diff>0.05);

    % Define GPS object
    reference_latitude = 40.86368573;
    reference_longitude = -77.83592832;
    reference_altitude = 344.189;
    gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

    % Use the class to convert LLA to ENU
    LLA_data = gps_object.ENU2WGSLLA(LiDAR_outer_edge);

    % Plot the ENU results
    figure(30000);

    % hold on


    geoplot(LLA_data(:,1),LLA_data(:,2),'mo','MarkerSize',10);
    hold on
    geoplot(LLA_data(:,1),LLA_data(:,2),'k.','MarkerSize',10);

    geoplot(LLA_data(aa,1),LLA_data(aa,2),'c.','MarkerSize',20);
    geoplot(LLA_data(aa,1),LLA_data(aa,2),'b.','MarkerSize',10);
    
    % geoplot(LLA_data(aa(end),1),LLA_data(aa(end),2),'c.','MarkerSize',20);
    % geoplot(LLA_data(aa(end),1),LLA_data(aa(end),2),'b.','MarkerSize',10);

    title('LLA Trace geometry')
        % xlabel('Latitude')
    % ylabel('Longitude')

    geobasemap satellite
    geotickformat -dd

end








