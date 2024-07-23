
%% script_test_geometry_updatedSurfaceAnalysis
% This script follows the updated surface analysis schema step by step

% Delete "fcn_geometry_classifyIntoGrids", "fcn_geometry_surfaceAnalysis"
% and "script_test_geometry_surface_analysis" after making this script into
% a function

% Revision History
% 2024_07_18 - Aneesh Batchu
% -- Wrote this code originally 

%% STEP 1: Load and study the data (Yet to be functionalized)

% This is done in "script_load_sample_LIDAR_data"
% This is done in "script_plot_sample_data"

% These scripts are in the main directory of Geom Class
% Just "run" the scripts. You don't need to run section by section

% concatenate_LiDAR_XYZ_points
% concatenate_scanLine_rings

%% STEP 2: Find the driven path (left and right side points) (Yet to be functionalized)

% This is done in "script_test_geometry_boundaryPointsDrivenPath" 

% This script can be found in "Functions" directory of "Geom Class" repo
% Just "run" the script. You don't need to run section by section

% boundary_points_driven_path
% boundary_points_driven_path_LLA

%% STEP 3 & STEP 4: Seperate the data into grids, and classify the grids as the grids with zero points and grids with more than zero points

% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
LiDAR_allPoints = concatenate_LiDAR_XYZ_points(:,1:3);

% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries in #D
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 

% plot all LiDAR points                     ----------- To-Do: ALEKS
%
% Change the function name into fcn_geometry_plotPointsinLLA
% - Add an option for chnaging the marker such as '.' (currently plots),
% '+', '*', 'o' etc
% BUG: This function does not work if Legend options are given as zero.
% Write some test cases in the script
% BUG: This function should work even if the legend_option is empty
% BUG: This function should work even if the Legend_name is empty
% Check this function carefully. Write all the test script clearly. It
% should cover all the cases. 
% 
% marker_size = 10;
% RGB_triplet = [0.8, 0.8, 0.8]; 
% legend_option = 1;
% legend_name = 'LiDAR Points';
% legend_position = [];
% plot_LiDAR_allPoints = LiDAR_allPoints; 
% [~] = fcn_geometry_plotGridCenters(LiDAR_allPoints,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% Define GPS object - This should be the input for the above function 
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
% Define GPS object
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert ENU to LLA
LIDAR_allPoints_LLA = gps_object.ENU2WGSLLA(LiDAR_allPoints);

% Currently plotting here without uisng any function 
fig_num_LLA = 30;
% Plot the ENU results
figure(fig_num_LLA);clf;

geoplot(LIDAR_allPoints_LLA(:,1),LIDAR_allPoints_LLA(:,2),'mo','MarkerSize',10);
hold on
geoplot(LIDAR_allPoints_LLA(:,1),LIDAR_allPoints_LLA(:,2),'k.','MarkerSize',10);
geoplot(boundary_points_driven_path_LLA(:,1),boundary_points_driven_path_LLA(:,2),'g.','MarkerSize',10);

title('LLA Trace geometry')

geobasemap satellite
geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default

fig_num_first_classification = 40; 
figure(fig_num_first_classification);clf

[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_geometry_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,fig_num_first_classification);

%% STEP 5: Find the driven path grids within the grids more than zero points

% -----------------------------NOTE------------------------------
% After finding the grids without anypoints, the grids are completely
% removed from the analysis. Only, grids with greater than zero points were
% analyzed from here. 
% -----------------------------NOTE------------------------------

% Plot all the grids greater than zero point density

% Figure number
fig_num_gridLines_greater_than_zero_point_density = 50;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')

% Pre-allocation: To find the total number of points in each grid
total_points_in_grids_greater_than_zero = zeros(length(grids_greater_than_zero_points),1);

% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11

total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

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

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(concatenate_scanLine_rings(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    
    % Find the total number of points in each grid
    total_points_in_grids_greater_than_zero(ith_grid) = length(points_in_domain);
    
    % Save the grid lines of all the grids greater than zero density in a
    % matrix
    length_gridlines = length(gridlines);
    gridlines_grids_greater_than_zero(1+(ith_grid-1)*length_gridlines:ith_grid*length_gridlines,:) = gridlines;
   
    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % plot the points in the grid
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.0d',ith_text);

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

%% STEP 5: Plot grid centers and boundary points in a same figure in ENU

% "inpolygon" is used to find the grids within the boundary points 
[in,on] = inpolygon(gridCenters_greater_than_zero_point_density(:,1),gridCenters_greater_than_zero_point_density(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = grids_greater_than_zero_points(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = find(in); 

% Total points in each grid in the driven path
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(grids_greater_than_zero_points); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2)];


fig_num_gridCenters_and_boundary_points_greater_than_zero = 51;
figure(fig_num_gridCenters_and_boundary_points_greater_than_zero); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_greater_than_zero_point_density(:,1), gridCenters_greater_than_zero_point_density(:,2), '.','MarkerSize',40,'Color',[0.8 0.8 0.8]);
plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.0d',ith_text);
    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 8);
end

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',20,'Color',[0 1 0]) % points strictly inside


%% STEP 6: Statistic 1

% Figure number of histogram
fig_histogram = 52; 
figure(fig_histogram); clf; 

% edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

% Create the histogram
% actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on'); 
actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,'Visible','on'); 
hold on 
% total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on'); 
total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,'Visible','on'); 

% Extract the counts for both histograms
counts1 = actual_driven_path_grids_hist.Values;
counts2 = total_grids_greater_than_zero_hist.Values;
binEdges = total_grids_greater_than_zero_hist.BinEdges;

% Calculate the overlapping counts
% overlapCounts = min(counts2, counts1);

% Find a ratio
point_density = sum(binEdges(1:2))/2; 

% Add labels and title 
xlabel('Points per grid');
ylabel('Frequency');
title('Histogram of points per grid');

%% STEP 6: Statistic 2 - Determine number of LiDAR scans in each grid

% This was done in STEP 5, however, it was done using a for loop. Need to
% do it without using a for loop.

% total_scan_lines_in_each_grid

%% STEP 7: Classify the grids with more than zero points into mapped and unmapped grids

% plotting
fig_num = 72; 
figure(fig_num);clf

% Classify the grids with more than zero points into mapped and unmapped grids
[original_grids_with_low_point_density, original_grids_with_required_point_density, original_grids_with_more_than_one_scan_line, original_grids_with_one_scan_line, ...
    original_mapped_grids, original_unmapped_grids, gridCenters_low_point_density, gridCenters_required_point_density, gridCenters_with_more_than_one_scan_line, ...
gridCenters_with_one_scan_line, gridCenters_mapped_grids, gridCenters_unmapped_grids, current_grids_with_low_point_density, current_grids_with_required_point_density, ...
current_grids_with_more_than_one_scan_line, current_grids_with_one_scan_line, current_mapped_grids, current_unmapped_grids]...
= fcn_geometry_GridsIntoMappedUnmapped(point_density, total_N_points_in_each_grid, total_scan_lines_in_each_grid_with_more_than_zero_points, ...
grids_greater_than_zero_points, gridCenters, fig_num); 

% Plot the grids with low point density and required density 
fig_num = 70; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with low point density';
legend_position = [];
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_low_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with required density';
legend_position = [];
plot_gridCenters_required_point_density = [gridCenters_required_point_density, zeros(length(gridCenters_required_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_required_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% Plot grids with one scan line and more than one scan line
fig_num = 71; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with more than one scan line';
legend_position = [];
plot_gridCenters_with_more_than_one_scan_line = [gridCenters_with_more_than_one_scan_line, zeros(length(gridCenters_with_more_than_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_with_more_than_one_scan_line,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% ENU - plotting mapped and unmapped
fig_num_ENU = 73; 
figure(fig_num_ENU)
clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of mapped and unmapped grids and boundary points')

plot(gridCenters_unmapped_grids(:,1), gridCenters_unmapped_grids(:,2), '.','MarkerSize',40,'Color',[0.8 0.8 0.8]);
plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_unmapped_grids(:,1))
    current_text = sprintf('%.0d',current_unmapped_grids(ith_text));
    % Place the text on the grid center
    text(gridCenters_unmapped_grids(ith_text,1), gridCenters_unmapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 8);
end

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8);
end


%% STEP 8: Statistic 3 - Standard deviation in Z

input_points = LiDAR_allPoints; 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_mapped_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

for ith_mapped_grid = 1:total_mapped_grids
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, ~, ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

fig_num_ENU_statistic_three = 80; 
figure(fig_num_ENU_statistic_three);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2) % points strictly inside


% Plot grid lines and standard deviation
fig_num = 81;
figure(fig_num)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')

% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(original_mapped_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(original_mapped_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];
    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_domain),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

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

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==original_mapped_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    % Plot the result
    % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % hold on
    % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
    gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;

end

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
standard_deviation_in_z_round = round(standard_deviation_in_z,3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_mapped_grids(:,1))
    current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end


% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(current_mapped_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
std_in_z_driven_path = standard_deviation_in_z(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
std_in_z_other_mapped_grids = standard_deviation_in_z(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 82;
figure(fig_num)

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(current_mapped_grids, standard_deviation_in_z,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, std_in_z_driven_path,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(current_mapped_grids(~driven_path_grid_indices_in_current_mapped_grids), std_in_z_other_mapped_grids,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_std_in_z_driven_path = mean(std_in_z_driven_path); 

% Find mean std in z of not driven path
max_std_in_z_not_driven_path = max(std_in_z_other_mapped_grids); 

% Std Threshold
% Instead of choosing 6, try to find a ratio
% ratio: mean_std_in_z_driven_path/mean_std_of_all_grids
std_threshold = mean_std_in_z_driven_path*6;

disp(std_threshold)


%% STEP 8: Statistic 4 - angle deviation

input_points = LiDAR_allPoints; 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_mapped_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


for ith_mapped_grid = 1:total_mapped_grids
    [~, ~, ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% STEP 2
% Comparing normal vector with vertical direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product);

fig_num_ENU_statistic_four = 83; 
figure(fig_num_ENU_statistic_four);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% Plot grid lines and standard deviation
fig_num = 84;
figure(fig_num)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Angle deviation of mapped grids')

% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(original_mapped_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(original_mapped_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];
    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_domain),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

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

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==original_mapped_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    % Plot the result
    % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % hold on
    % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
    gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;

end

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
angle_btw_unit_normals_and_vertical_round = round((angle_btw_unit_normals_and_vertical*180/pi),3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_mapped_grids(:,1))
    current_text = sprintf('%.3f',angle_btw_unit_normals_and_vertical_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end

% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(current_mapped_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
angle_btw_unit_normals_and_vertical_driven_path = angle_btw_unit_normals_and_vertical(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
angle_btw_unit_normals_and_vertical_other_mapped_grids = angle_btw_unit_normals_and_vertical(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 85;
figure(fig_num)

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(current_mapped_grids, angle_btw_unit_normals_and_vertical*180/pi,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, angle_btw_unit_normals_and_vertical_driven_path*180/pi,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(current_mapped_grids(~driven_path_grid_indices_in_current_mapped_grids), angle_btw_unit_normals_and_vertical_other_mapped_grids*180/pi,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path); 

% Find mean std in z of not driven path
max_angle_btw_unit_normals_and_vertical_not_driven_path = max(angle_btw_unit_normals_and_vertical_other_mapped_grids); 

% Theta threshold
% RATIO: Find a ratio
theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 0.04;

disp(theta_threshold*(180/pi))


%% STEP 9: 3rd Classification

input_points = LiDAR_allPoints; 
% std_threshold = 0.05; 
% theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters

fig_num = 6000; 
figure(fig_num);clf


% Classify mapped grids into drivable and drivable
[standard_deviation_in_z, angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids, current_non_drivable_grid_numbers_in_mapped_grids, ...
    gridCenters_drivable_grids,gridCenters_non_drivable_grids, concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_geometry_classifyGridsAsDrivable(gridIndices_cell_array, original_mapped_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num); 


%% STEP 10: Find the boundary points of drivable and non-drivable grids

% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_mapped_unmapped = 767787; 

XYZ_matrix_mapped_grids = [gridCenters_mapped_grids(:,1:2) ones(length(gridCenters_mapped_grids(:,1)),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_unmapped_grids(:,1:2) zeros(length(gridCenters_unmapped_grids(:,1)),1)]; 

XYZ_matrix_mapped_unmapped_gridcenters = [XYZ_matrix_mapped_grids; XYZ_matrix_unmapped_grids]; 

% Find the unique elements
XYZ_matrix = unique(XYZ_matrix_mapped_unmapped_gridcenters,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4);

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

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_mapped_unmapped);
clf;
boundary_points_mapped_unmapped = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_mapped_unmapped);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;

% XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)];
XYZ_matrix = concatenate_gridCenters_drivable_non_drivable_grids;

XYZ_matrix = unique(XYZ_matrix,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

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

figure(fig_num_drivable_non_drivable)
clf;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_drivable_non_drivable);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 

fig_num_bd_pts_ENU = 1000; 

[members, id_x] = ismember(boundary_points,boundary_points_mapped_unmapped,'rows'); 

not_boundary_points = boundary_points(members,:);

true_boundary_points = boundary_points(members==0,:);

figure(fig_num_bd_pts_ENU)
clf;
% figure(10001)
hold on
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40)
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30)

fig_num = 6000;
% plot computed boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotGridCenters(plot_true_boundary_points,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot computed boundary points
marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Computed boundary points';
legend_position = [];

[~] = fcn_geometry_plotGridCenters(plot_true_boundary_points,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);


% plot driven path
marker_size = 30;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];

plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_driven_path,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% Find the nearest boundary points in ENU - DrB

% Instructions
% Run this script, after running the scripts in STEP 1 and STEP 2. 


% Write the grid number at the grid center for reference. 

fig_num = 3332; 
figure(fig_num);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of drivable and non-drivable grids in ENU')

% gridCenters_drivable_grids
% 
% gridCenters_non_drivable_grids

% plot(gridCenters_mapped_grids(:,1), )

% plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);
p1 = plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',45,'Color',[0.4660 0.6740 0.1880],'DisplayName','Drivable grids');
p2 = plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',45,'Color',[0.6350 0.0780 0.1840], 'DisplayName','Non-drivable grids');

for ith_text = 1:length(original_mapped_grids(:,1))
    current_text = sprintf('%.0d',original_mapped_grids(ith_text));

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1, 1, 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% plot true boundary points
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
p3 = plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');

% % plot the grids in the driven path
p4 = plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2, 'DisplayName','Driven path grids'); % points strictly inside
% % plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 0 0], 'LineWidth',0.5) % points strictly inside

% ----------------------- TRAVERSING TOP TO BOTTOM -------------------------
% Create a combined grid center matrix of mapped grid centers and boundary
% points
combined_gridCenters = [gridCenters_mapped_grids; true_boundary_points];

% Indices of the sorted the combined Grid Centers in y direction in
% descending order
[~, sorted_combined_gridCenters_indices_in_y_direction] = sort(combined_gridCenters(:,2),'descend'); 

% Sorted combined grid centers in y direction in descending order
sorted_combined_gridCenters_in_y_direction = combined_gridCenters(sorted_combined_gridCenters_indices_in_y_direction, :); 

% Round the sorted combined gridcenters to the fourth decimal
sorted_combined_gridCenters_in_y_direction = round(sorted_combined_gridCenters_in_y_direction,4); 

% Indices of the sorted grid centers of the driven path in x direction
[~, sorted_driven_path_gridCenters_indices_in_x_direction] = sort(gridCenters_driven_path(:,1));

% Sorted driven path grid centers in x direction
sorted_driven_path_gridCenters_in_x_direction = gridCenters_driven_path(sorted_driven_path_gridCenters_indices_in_x_direction, :); 

% Round the sorted driven path gridcenters to the fourth decimal
sorted_driven_path_gridCenters_in_x_direction = round(sorted_driven_path_gridCenters_in_x_direction,4); 

% Find unique x coordinates of driven path grid centers
[unique_x_coord_driven_path, unique_x_coor_driven_path_index] = unique(sorted_driven_path_gridCenters_in_x_direction(:,1)); 

% Unique grid centers of driven path sorted in x direction
driven_path_gridCenters_sorted_in_x_direction = sorted_driven_path_gridCenters_in_x_direction(unique_x_coor_driven_path_index,:);

% Declare the top and bottom boundary points matrices
top_boundary_points = []; 
bottom_boundary_points = []; 

% This loop traverse the points from top to bottom 
for ith_grid_center = 1:length(unique_x_coord_driven_path)

    % Store the indices of grid centers whose x coordinates match with the
    % x coordinates of driven path
    matching_gridIndices = find(sorted_combined_gridCenters_in_y_direction(:,1) == unique_x_coord_driven_path(ith_grid_center)); 
    
    % Get the matched grid centers using the stored indices
    matched_gridCenters = sorted_combined_gridCenters_in_y_direction(matching_gridIndices,:); 

    % Find the difference of y coordinate of driven path grid center
    % sorted in x direction from y coordinates of matched grid centers
    y_coordinate_difference = matched_gridCenters(:,2) - driven_path_gridCenters_sorted_in_x_direction(ith_grid_center, 2);

    % Divide the above output with grid size and multiply with 2. The odd
    % number in the output indicates the grid center
    scaled_y_coordinate_difference = (y_coordinate_difference/grid_size)*2; 
    
    % Total number of top grid centers + bpoundary points
    total_top_gridCenters = sum(scaled_y_coordinate_difference > 0);
    
    % Traverse top grids from the driven path grid center
    for top_grid = total_top_gridCenters:-1:1
        % Find the odd number in the top preliminary output. That grid
        % index indicate the grid center of the top boundary point
        if (mod(scaled_y_coordinate_difference(top_grid),2) == 1)
            % The top boundary point index from the top
            top_boundary_point_index = top_grid;
            % Top boundary points
            top_boundary_points = [top_boundary_points; matched_gridCenters(top_boundary_point_index,:)]; %#ok<AGROW>
            break
        end
        top_boundary_point_index = NaN;  
    end

    % Total number of bottom grids
    total_bottom_gridCenters = sum(scaled_y_coordinate_difference < 0);

    % Traverse bottom grids from the driven path grid center
    for bottom_grid = 1:total_bottom_gridCenters
        % Find the odd number in the top preliminary output. That grid
        % index indicate the grid center of the bottom boundary point
        bottom_index_offset = total_top_gridCenters + 1; 
        if (mod(scaled_y_coordinate_difference(bottom_index_offset+bottom_grid),2) == 1)
            % The bottom boundary point index from the top
            bottom_boundary_point_index = bottom_index_offset + bottom_grid;
            % bottom boundary points
            bottom_boundary_points = [bottom_boundary_points; matched_gridCenters(bottom_boundary_point_index,:)]; %#ok<AGROW>
            break
        end
        bottom_boundary_point_index = NaN; 
    end

end

% plot the top-bottom boundary points
p5 = plot(top_boundary_points(:,1), top_boundary_points(:,2), 'o','MarkerSize',20,'Color',[1 0 0], 'LineWidth',2, 'DisplayName','Nearest boundary points');
p6 = plot(bottom_boundary_points(:,1), bottom_boundary_points(:,2), 'o','MarkerSize',20,'Color',[1 0 0], 'LineWidth',2, 'DisplayName','Nearest boundary points'); 


% ----------------------- TRAVERSING LEFT TO RIGHT -------------------------

% Indices of the sorted the combined Grid Centers in x direction in
% ascending order
[~, sorted_combined_gridCenters_indices_in_x_direction] = sort(combined_gridCenters(:,1)); 

% Sorted combined grid centers in x direction in descending order
sorted_combined_gridCenters_in_x_direction = combined_gridCenters(sorted_combined_gridCenters_indices_in_x_direction, :); 

% Round the sorted combined gridcenters to the fourth decimal
sorted_combined_gridCenters_in_x_direction = round(sorted_combined_gridCenters_in_x_direction,4); 

% Indices of the sorted grid centers of the driven path in y direction
[~, sorted_driven_path_gridCenters_indices_in_y_direction] = sort(gridCenters_driven_path(:,2));

% Sorted driven path grid centers in y direction
sorted_driven_path_gridCenters_in_y_direction = gridCenters_driven_path(sorted_driven_path_gridCenters_indices_in_y_direction, :); 

% Round the sorted driven path gridcenters to the fourth decimal
sorted_driven_path_gridCenters_in_y_direction = round(sorted_driven_path_gridCenters_in_y_direction,4); 

% Find unique y coordinates of driven path grid centers
[unique_y_coord_driven_path, unique_y_coor_driven_path_index] = unique(sorted_driven_path_gridCenters_in_y_direction(:,2)); 

% Unique grid centers of driven path sorted in x direction
driven_path_gridCenters_sorted_in_y_direction = sorted_driven_path_gridCenters_in_y_direction(unique_y_coor_driven_path_index,:);

% Declare the top and bottom boundary points matrices
left_boundary_points = []; 
right_boundary_points = []; 

% This loop traverse the points from left to right 
for ith_grid_center = 1:length(unique_y_coord_driven_path)

    % Store the indices of grid centers whose y coordinates match with the
    % y coordinates of driven path
    matching_gridIndices = find(sorted_combined_gridCenters_in_x_direction(:,2) == unique_y_coord_driven_path(ith_grid_center)); 
    
    % Get the matched grid centers using the stored indices
    matched_gridCenters = sorted_combined_gridCenters_in_x_direction(matching_gridIndices,:); 

    % Find the difference of x coordinate of driven path grid center
    % sorted in y direction from x coordinates of matched grid centers
    x_coordinate_difference = matched_gridCenters(:,1) - driven_path_gridCenters_sorted_in_y_direction(ith_grid_center, 1);

    % Divide the above output with grid size and multiply with 2. The odd
    % number in the output indicates the grid center
    scaled_x_coordinate_difference = (x_coordinate_difference/grid_size)*2; 
    
    % Total number of top grid centers + bpoundary points
    total_left_gridCenters = sum(scaled_x_coordinate_difference < 0);
    
    % Traverse top grids from the driven path grid center
    for left_grid = total_left_gridCenters:-1:1
        % Find the odd number in the top preliminary output. That grid
        % index indicate the grid center of the top boundary point
        if (mod(scaled_x_coordinate_difference(left_grid),2) == 1)
            % The top boundary point index from the top
            left_boundary_point_index = left_grid;
            % Top boundary points
            left_boundary_points = [left_boundary_points; matched_gridCenters(left_boundary_point_index,:)]; %#ok<AGROW>
            break
        end
        left_boundary_point_index = NaN; 
    end

    % Total number of bottom grids
    total_right_gridCenters = sum(scaled_x_coordinate_difference > 0);

    % Traverse bottom grids from the driven path grid center
    for right_grid = 1:total_right_gridCenters
        % Find the odd number in the top preliminary output. That grid
        % index indicate the grid center of the bottom boundary point
        right_index_offset = total_left_gridCenters + 1; 
        if (mod(scaled_x_coordinate_difference(right_index_offset+right_grid),2) == 1)
            % The bottom boundary point index from the top
            right_boundary_point_index = right_index_offset + right_grid;
            % bottom boundary points
            right_boundary_points = [right_boundary_points; matched_gridCenters(right_boundary_point_index,:)]; %#ok<AGROW>
            break
        end
        right_boundary_point_index = NaN; 
    end

end

% plot the left-right boundary points
if ~isempty(left_boundary_points)
    plot(left_boundary_points(:,1), left_boundary_points(:,2), 'o','MarkerSize',20,'Color',[1 0 0], 'LineWidth',2)
end
p7 = plot(right_boundary_points(:,1), right_boundary_points(:,2), 'o','MarkerSize',20,'Color',[1 0 0], 'LineWidth',2,'DisplayName','Nearest boundary points'); 




legend([p1 p2 p3 p4 p5])






