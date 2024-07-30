
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
% This is done in "script_plot_sample_data" - Functionalize this first (Aleks)

% These scripts are in the main directory of Geom Class
% Just "run" the scripts. You don't need to run section by section

% concatenate_LiDAR_XYZ_points
% concatenate_scanLine_rings

%% STEP 2: Find the driven path (left and right side points) (Yet to be functionalized)

% This is done in "script_test_geometry_boundaryPointsDrivenPath" - Steven 

% This script can be found in "Functions" directory of "Geom Class" repo
% Just "run" the script. You don't need to run section by section

% boundary_points_driven_path
% boundary_points_driven_path_LLA


%% STEP 3 & STEP 4: Seperate the data into grids, and classify the grids as the grids with zero points and grids with more than zero points

% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
% LiDAR_allPoints = concatenate_LiDAR_XYZ_points(:,1:3);

LiDAR_allPoints = [concatenate_LiDAR_XYZ_points, concatenate_scanLine_rings];

% remove NANs
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);

% scanLine_rings without NaNs
LIDAR_scanLines = LiDAR_allPoints(:,4:5); 


% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 0.3; %0.8;%1;%1.25; %1.26

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
LIDAR_allPoints_LLA = gps_object.ENU2WGSLLA(LiDAR_allPoints(:,1:3));

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

concatenate_gridPoints_scanLines_rings = [];

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
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
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];
    % gridPoints_scanLines_rings_to_add = [gridPoints_scanLines_rings_to_add; naa nan nan nan]; %#ok<AGROW>

    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>


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

%% Testing: Delete this afterwards


% Figure number
fig_num_gridLines_greater_than_zero_point_density = 545;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')


% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];
orthogonal_dist_each_grid = []; 

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
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
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];

    % Sort gridPoints_scanLines_rings_to_add based on scan line
    [~,sorted_gridPoints_scanLines_rings] = sort(gridPoints_scanLines_rings_to_add(:,4));

    % Sorted gridPoints_scanLines_rings_to_add matrix
    gridPoints_scanLines_rings_to_add = gridPoints_scanLines_rings_to_add(sorted_gridPoints_scanLines_rings,:);

    % Indices first scan line of the matrix as a seperate matrix
    indices_gridPoints_scanLines_first_scan = find(gridPoints_scanLines_rings_to_add(:,4) == gridPoints_scanLines_rings_to_add(1,4)); 

    % Seperate the first scan line of the matrix as a seperate matrix
    gridPoints_scanLines_first_scan = gridPoints_scanLines_rings_to_add(indices_gridPoints_scanLines_first_scan,:);

    %
    if length(gridPoints_scanLines_first_scan(:,1)) > 1 && (gridPoints_scanLines_first_scan(1,5) == gridPoints_scanLines_first_scan(2,5))
        change_in_vector = gridPoints_scanLines_first_scan(2,2:3) - gridPoints_scanLines_first_scan(1,2:3);
        unit_change_in_vector = fcn_geometry_calcUnitVector(change_in_vector);
        orth_unit_change_in_vector = unit_change_in_vector*[0 1; -1 0];
        
        % The remaining number of grids
        remaining_grids = length(indices_gridPoints_scanLines_first_scan)+1:length(gridPoints_scanLines_rings_to_add(:,1)); 
        
        % 
        vector_from_base_point_first_scan_to_points_in_otherScans_rings = gridPoints_scanLines_rings_to_add(remaining_grids,2:3) - ...
            gridPoints_scanLines_first_scan(1,2:3).*(ones(length(remaining_grids),2));

        % Unit orthogonal vector
        repeated_orth_unit_change_in_vector = orth_unit_change_in_vector.*(ones(length(remaining_grids),2));

        % Calculate the transverse distance
        transverse_dist_grid_points_other_scanLines = sum(vector_from_base_point_first_scan_to_points_in_otherScans_rings.*repeated_orth_unit_change_in_vector,2);

        % Mean of absolute values of transverse distances
        mean_dist = mean(abs(transverse_dist_grid_points_other_scanLines));

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; mean_dist]; %#ok<AGROW>
    else 

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; nan]; %#ok<AGROW>

    end

   

    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>

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
    % current_text = sprintf('%.0d',current_mapped_grids(ith_text));
    current_text = sprintf('%.0d',(ith_text));

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

fig_num = 5837; 
figure(fig_num); clf;
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
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

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.3f',orthogonal_dist_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

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

grid_corners = grid_AABBs(grids_greater_than_zero_points(in),1:4); 

% Extract the coordinates
x_low = grid_corners(:,1); 
x_high = grid_corners(:,2);
y_low = grid_corners(:,3);
y_high =  grid_corners(:,4); 

% Rearrange grid corners
rearranged_grid_corners = [x_low y_low; x_high y_low; x_low y_high; x_high y_high]; 


% "inpolygon" is used to find the grid corners within the boundary points 
[in_bb,on_bb] = inpolygon(rearranged_grid_corners(:,1),rearranged_grid_corners(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));


inliers_corners = rearranged_grid_corners(in_bb,:); 


fig_num_gridCenters_and_boundary_points_greater_than_zero = 51;
figure(fig_num_gridCenters_and_boundary_points_greater_than_zero); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_greater_than_zero_point_density(:,1), gridCenters_greater_than_zero_point_density(:,2), '.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 
plot(rearranged_grid_corners(in_bb,1), rearranged_grid_corners(in_bb,2), '.', 'MarkerSize',10, 'Color',[0 1 0]); 

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5) % points strictly inside
aaa = grid_AABBs(grids_greater_than_zero_points(in),1:4);

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.0d',ith_text);
    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 6);
end


%% STEP 6: Statistic 1

% Figure number of histogram
fig_histogram = 52; 
figure(fig_histogram); clf; 
% Add labels and title 
hold on
grid on
xlabel('Points per grid');
ylabel('Frequency');
title('Statistic 1: Determining suitable point density');


% edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

% Create the histogram
% actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on'); 

% total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on'); 
total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,20,'Visible','on'); 
actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,10,'Visible','on'); 

% Extract the counts for both histograms
counts1 = total_grids_greater_than_zero_hist.Values;
[~,index_max_counts2] = max(counts1);
counts2 = actual_driven_path_grids_hist.Values;

binEdges = total_grids_greater_than_zero_hist.BinEdges;

% Calculate the overlapping counts
% % overlapCounts = min(counts2, counts1);

% Find a ratio
 % point_density = sum(binEdges(index_max_counts2:(index_max_counts2+1)))/2; 

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 1.5*(std(total_points_in_each_grid_in_the_driven_path)));

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path));


% Minimum number of points required 
% point_density =
% mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);

plot(point_density,0, 'k.', 'MarkerSize',20)
current_text = sprintf('point density = %.2d',point_density);
text(650, 200,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');


%% STEP 6: Statistic 2 - Determine number of LiDAR scans in each grid

% This was also done in STEP 5, however, it was done using a for loop. Need to
% do it without using a for loop.

total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
     %Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);

     %Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
  
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORK     S WITHOUT A FOR LOOP
    total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
end


% total_scan_lines_in_each_grid
%wrote by Jiabao Zhao

%total_scan_lines_in_each_grid_with_more_than_zero_points = [];
%rows_in_domain = ismember(gridIndices, grids_greater_than_zero_points);
%scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
%total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid];


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
marker_type = [];
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_low_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with required density';
legend_position = [];
marker_type = [];
plot_gridCenters_required_point_density = [gridCenters_required_point_density, zeros(length(gridCenters_required_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_required_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% Plot grids with one scan line and more than one scan line
fig_num = 71; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with more than one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_more_than_one_scan_line = [gridCenters_with_more_than_one_scan_line, zeros(length(gridCenters_with_more_than_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_more_than_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


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


%% STEP 7.25: Find the orthogonal distances 

% Figure number
fig_num_gridLines_greater_than_zero_point_density = 50;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')


% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(original_mapped_grids),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];
orthogonal_dist_each_grid = []; 

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_grid),1:4);

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
    rows_in_domain = gridIndices==original_mapped_grids(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];

    % Sort gridPoints_scanLines_rings_to_add based on scan line
    [~,sorted_gridPoints_scanLines_rings] = sort(gridPoints_scanLines_rings_to_add(:,4));

    % Sorted gridPoints_scanLines_rings_to_add matrix
    gridPoints_scanLines_rings_to_add = gridPoints_scanLines_rings_to_add(sorted_gridPoints_scanLines_rings,:);

    % Count occurrences of each unique number in scan lines
    [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_rings_to_add(:,4));
    counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

    % Create the new array with the counts
    count_array = counts;

    % Index of the scan line with more than one occurence
    index_of_scanLines = find(count_array>1, 1);

    % Indices first scan line of the matrix as a seperate matrix
    indices_gridPoints_scanLines_first_scan = find(gridPoints_scanLines_rings_to_add(:,4) == gridPoints_scanLines_rings_to_add(index_of_scanLines,4));

    % Seperate the scan line with more than one occurence of the matrix as a seperate matrix
    gridPoints_scanLines_first_scan = gridPoints_scanLines_rings_to_add(indices_gridPoints_scanLines_first_scan,:);

    % Count occurrences of each unique number in rings
    [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_first_scan(:,5));
    counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

    % Create the new array with the counts
    count_array = counts;

    % Index of the scan line with more than one occurence
    index_of_rings = find(count_array>1, 1);

    % if length(gridPoints_scanLines_first_scan(:,1)) == 1
    % 
    %     indices_gridPoints_scanLines_first_scan

    if length(gridPoints_scanLines_first_scan(:,1)) > 1 && (gridPoints_scanLines_first_scan(index_of_rings,5) == gridPoints_scanLines_first_scan(index_of_rings+1,5))
        change_in_vector = gridPoints_scanLines_first_scan(2,2:3) - gridPoints_scanLines_first_scan(1,2:3);
        unit_change_in_vector = fcn_geometry_calcUnitVector(change_in_vector);
        orth_unit_change_in_vector = unit_change_in_vector*[0 1; -1 0];
        
        % The remaining number of grids
        remaining_grids = length(indices_gridPoints_scanLines_first_scan)+1:length(gridPoints_scanLines_rings_to_add(:,1)); 
        
        % 
        vector_from_base_point_first_scan_to_points_in_otherScans_rings = gridPoints_scanLines_rings_to_add(remaining_grids,2:3) - ...
            gridPoints_scanLines_first_scan(1,2:3).*(ones(length(remaining_grids),2));

        % Unit orthogonal vector
        repeated_orth_unit_change_in_vector = orth_unit_change_in_vector.*(ones(length(remaining_grids),2));

        % Calculate the transverse distance
        transverse_dist_grid_points_other_scanLines = sum(vector_from_base_point_first_scan_to_points_in_otherScans_rings.*repeated_orth_unit_change_in_vector,2);

        % Mean of absolute values of transverse distances
        mean_dist = mean(abs(transverse_dist_grid_points_other_scanLines));

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; mean_dist]; %#ok<AGROW>
    else 

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; nan]; %#ok<AGROW>

    end

    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>

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

for ith_text = 1:length(original_mapped_grids(:,1))
    % current_text = sprintf('%.0d',current_mapped_grids(ith_text));
    current_text = sprintf('%.0d',(ith_text));

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

fig_num = 5837; 
figure(fig_num); clf;
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(original_mapped_grids)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_grid),1:4);

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

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.3f',orthogonal_dist_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

%% Updated mapped grids

threshold_transverse_dist = 0.03;

updated_original_mapped_grids = original_mapped_grids(orthogonal_dist_each_grid>threshold_transverse_dist); 

updated_original_unmapped_grids = sort([original_unmapped_grids; original_mapped_grids(orthogonal_dist_each_grid<=threshold_transverse_dist)]); 

gridCenters_updated_original_mapped_grids = gridCenters(updated_original_mapped_grids,1:2); 

gridCenters_updated_original_unmapped_grids = gridCenters(updated_original_unmapped_grids,1:2); 

% Plot the grids with low point density and required density 
fig_num = 8765; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Updated mapped grids';
legend_position = [];
marker_type = [];
plot_gridCenters_updated_original_mapped_grids = [gridCenters_updated_original_mapped_grids, zeros(length(gridCenters_updated_original_mapped_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_mapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Updated unmapped grids';
legend_position = [];
marker_type = [];
plot_gridCenters_updated_original_unmapped_grids = [gridCenters_updated_original_unmapped_grids, zeros(length(gridCenters_updated_original_unmapped_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% Find mapped grids current number

% ENU - plotting mapped and unmapped
fig_num_ENU = 74; 
figure(fig_num_ENU)
clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of mapped grids')

plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside



for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end

%% STEP 7.5: Find driven path of the vehicle in mapped grids

% "inpolygon" is used to find the grids within the boundary points 
[in,on] = inpolygon(gridCenters_mapped_grids(:,1),gridCenters_mapped_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = original_mapped_grids(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = current_mapped_grids(in); %find(in); 

% Total points in each grid in the driven path
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(current_mapped_grids); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_mapped_grids(in,1),gridCenters_mapped_grids(in,2)];


fig_num_gridCenters_and_boundary_points_greater_than_zero = 51;
figure(fig_num_gridCenters_and_boundary_points_greater_than_zero); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5) % points strictly inside

for ith_text = 1:length(gridCenters_mapped_grids(:,1))
    current_text = sprintf('%.0d',ith_text);
    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 6);
end


%% STEP 8: Statistic 3 - Standard deviation in Z

input_points =  LiDAR_allPoints(:,1:3);

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_mapped_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

for ith_mapped_grid = 1:total_mapped_grids
    % points = input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:);
    % points = points(~isnan(points(:,1)),:);
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

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside


for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end



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
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

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
mean_std_in_z_not_driven_path = mean(std_in_z_other_mapped_grids(~isnan(std_in_z_other_mapped_grids))); 

% Find max std in z of not driven path
max_std_in_z_not_driven_path = max(std_in_z_other_mapped_grids); 

% Std Threshold
% Instead of choosing 6, try to find a ratio
% ratio: mean_std_in_z_driven_path/mean_std_of_all_grids
% std_threshold = mean_std_in_z_driven_path*6;
% 
% disp(std_threshold)

%% Histogram of standard deviation

figure(123);clf;
hold on
grid on
xlabel('Standard deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable standard deviation in z')
% histogram(standard_deviation_in_z)
histogram(standard_deviation_in_z,100)
histogram(std_in_z_driven_path,5)

% std_threshold = mean_std_in_z_driven_path + 6*std(std_in_z_driven_path); 
% std_threshold = mean_std_in_z_driven_path + 3*std(std_in_z_driven_path); 
std_threshold = 0.06; 
% plot(std_threshold,max(std_in_z_driven_path), 'k.', 'MarkerSize',20)
disp('mean of std_threshold of driven path')
disp(mean_std_in_z_driven_path)
disp('chosen std_threshold')
disp(std_threshold)
% std_threshold = 0.06; 
plot(std_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('std threshold = %.4f',std_threshold);
text(0.4, 80,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');


%% STEP 8: Statistic 4 - angle deviation

input_points = LiDAR_allPoints(:,1:3); 

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
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

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
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

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
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 0.04;
% 
% disp(theta_threshold*(180/pi))

%% Histogram of angle deviation

figure(1223);clf;
hold on
grid on
xlabel('Angle deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable angle deviation')
% histogram(standard_deviation_in_z)
histogram(angle_btw_unit_normals_and_vertical,100)
histogram(angle_btw_unit_normals_and_vertical_driven_path,5)

% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.7*std(angle_btw_unit_normals_and_vertical_driven_path); 
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.3*std(angle_btw_unit_normals_and_vertical_driven_path); 

% plot(theta_threshold,max(angle_btw_unit_normals_and_vertical), 'k.', 'MarkerSize',20)
% theta_threshold = 0.1507; 
theta_threshold = 0.3; 
disp('mean of angle_deviation_driven_path')
disp(mean_angle_btw_unit_normals_and_vertical_driven_path*(180/pi))
disp('theta threshold')
disp(theta_threshold*(180/pi))

plot(theta_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('theta threshold = %.1f',theta_threshold*(180/pi));
text(0.5, 30,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');


%% STEP 9: 3rd Classification

input_points = LiDAR_allPoints(:,1:3); 
% std_threshold = 0.05; 
% theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters

fig_num = 6000; 
figure(fig_num);clf

% theta_threshold = [];
% std_threshold = [];

% Classify mapped grids into drivable and drivable
[standard_deviation_in_z, angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids, current_non_drivable_grid_numbers_in_mapped_grids, ...
    gridCenters_drivable_grids,gridCenters_non_drivable_grids, concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_geometry_classifyGridsAsDrivable(gridIndices_cell_array, original_mapped_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

std_threshold_failed_gridCenters = gridCenters_mapped_grids(standard_deviation_in_z>std_threshold,:); 

theta_threshold_failed_gridCenters = gridCenters_mapped_grids(angle_btw_unit_normals_and_vertical>theta_threshold,:); 

% plot(std_threshold_failed_gridCenters(:,1), std_threshold_failed_gridCenters(:,2), 'o','MarkerSize',20,'Color',[0 1 1], 'LineWidth',2) 

% fig_num = 98989; 
% figure(fig_num); clf; 

marker_size = 10;
RGB_triplet = [1, 1, 0]; 
legend_option = 1;
legend_name = 'Std threshold failed grids';
marker_type = 'o';
legend_position = [];
plot_std_threshold_failed_gridCenters = [std_threshold_failed_gridCenters, zeros(length(std_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_std_threshold_failed_grids] = fcn_geometry_plotPointsinLLA(plot_std_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% fcn_geometry_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

% geoplot(LLA_data_std_threshold_failed_grids(:,1), LLA_data_std_threshold_failed_grids(:,2), 'o','MarkerSize',10,'Color',[1 1 0], 'LineWidth',2) 

% plot grid centers
marker_size = 15;
RGB_triplet = [1, 0, 1]; 
legend_option = 1;
legend_name = 'Theta threshold failed grids';
marker_type = 'o';
legend_position = [];
plot_theta_threshold_failed_gridCenters = [theta_threshold_failed_gridCenters, zeros(length(theta_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_theta_threshold_failed_grids] = fcn_geometry_plotPointsinLLA(plot_theta_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% geoplot(LLA_data_theta_threshold_failed_grids(:,1), LLA_data_theta_threshold_failed_grids(:,2), 'o','MarkerSize',15,'Color',[1 0 1], 'LineWidth',2) 


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
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot computed boundary points
marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 30;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% Find the nearest boundary points in ENU - updated

% grid_width = ceil(max(true_boundary_points(:,1)) - min(true_boundary_points(:,1))); 
% grid_height = ceil(max(true_boundary_points(:,2)) - min(true_boundary_points(:,2))); 
% plot computed boundary points

fig_num = 6020;
% plot computed boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Computed boundary points';
legend_position = [];

[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);



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

% for ith_text = 1:length(original_mapped_grids(:,1))
%     current_text = sprintf('%.0d',original_mapped_grids(ith_text));
% 
%     % Place the text on the grid center
%     text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1, 1, 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
% end

% plot true boundary points
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
p3 = plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');

% % plot the grids in the driven path
p4 = plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2, 'DisplayName','Driven path grids'); % points strictly inside
% % plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 0 0], 'LineWidth',0.5) % points strictly inside

% %% Remove boundary points on the left
% 
% % scanLineNumber_start = 1400; 
% % scanLineNumber_end = 1410; 
% 
% boundaryLineNumber_start = scanLineNumber_start - 8;
% boundaryLineNumber_end = scanLineNumber_end - 6; 
% 
% VehiclePose_current = VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:2); 
% % Calculate the vehicle orientation
% vehicle_change_in_pose_XY = diff(VehiclePose_current(:,1:2));
% 
% % Repeat the last value again, since diff removes one row. We want the same
% % number of vectors as the number of points, and diff removed one point.
% vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];
% 
% % Convert these to unit vectors
% unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);
% 
% % Find orthogonal vetors by rotating by 90 degrees in the CCW direction
% unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];
% 
% % Calculate the vectors
% vector_from_vehicle_pose_to_boundary_points = true_boundary_points - VehiclePose_current(28,:).*(ones(length(true_boundary_points(:,1)),2)); 
% 
% % Unit orthogonal vector
% repeated_unit_ortho_vehicle_vectors_XYZ = unit_ortho_vehicle_vectors_XY(28,:).*(ones(length(true_boundary_points(:,1)),2)); 
% 
% 
% % Calculate the transverse distance
% transverse_dist_boundary_points = sum(vector_from_vehicle_pose_to_boundary_points.*repeated_unit_ortho_vehicle_vectors_XYZ,2);
% 
% % Transverse distances of the right boundaries
% transverse_dist_right_boundary_points = transverse_dist_boundary_points(transverse_dist_boundary_points<0,:);
% 
% % Boundary points on the right
% boundary_points_right = true_boundary_points(transverse_dist_boundary_points<0,:);
% 
% % nearest boundary points to the right
% nearest_boundary_points_right = boundary_points_right(abs(transverse_dist_right_boundary_points)<2.8,:);
% 
% figure(3543);clf; 
% hold on
% axis on
% grid on 
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Right Points')
% 
% plot(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1), VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),'.','Color',[0 0 0],'MarkerSize',30) 
% 
% plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
% plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% plot(nearest_boundary_points_right(:,1), nearest_boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% 
% quiver(...
%     VehiclePose_current(:,1),VehiclePose_current(:,2),...
%     unit_vehicle_change_in_pose_XY(:,1),unit_vehicle_change_in_pose_XY(:,2),'-','LineWidth',3,'Color',[0 1 0]);
% 
% quiver(...
%     VehiclePose_current(:,1),VehiclePose_current(:,2),...
%     unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);
% % 
% % % quiver3(...
% % %     VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
% % %     unit_ortho_vehicle_vectors_XY(scanLineNumber_start,1),unit_ortho_vehicle_vectors_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
% % 
% % 
% % 
% %% plot boundary points to the right
% 
% 
% % grid_width = ceil(max(true_boundary_points(:,1)) - min(true_boundary_points(:,1))); 
% % grid_height = ceil(max(true_boundary_points(:,2)) - min(true_boundary_points(:,2))); 
% % plot computed boundary points
% 
% fig_num = 70067;clf;
% % plot computed boundary points
% marker_size = 30;
% RGB_triplet = [0 1 1]; 
% legend_option = 1;
% legend_name = 'Computed boundary points';
% legend_position = [];
% plot_nearest_boundary_points_right = [nearest_boundary_points_right, zeros(length(nearest_boundary_points_right),1)];
% [~] = fcn_geometry_plotGridCenters(plot_nearest_boundary_points_right,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);
% 
% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% 
% [~] = fcn_geometry_plotGridCenters(plot_nearest_boundary_points_right,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);
% 
% 
%%

oriGIN = [0 0 0]; 
unit_vertical = [0 0 1]; 
input_points = LiDAR_allPoints;
total_non_drivable_grids = length(original_non_drivable_grids);
for ith_mapped_grid = 7
    figure(ith_mapped_grid*10000); clf;

    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors, base_point, ~] =...
        fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),ith_mapped_grid*10000);
    xlabel('X[m]')
    ylabel('Y[m]')
    zlabel('Z[m]')
     title('Plane fit of a mapped grid')
    % title(sprintf('Plane fit of  %dth mapped grid',ith_mapped_grid))
    figure(ith_mapped_grid*11100)
    hold on;
    grid on;
    xlabel('X[m]')
    ylabel('Y[m]')
    zlabel('Z[m]')
    title(sprintf('angle btw %dth mapped grid (Non-drivable)',ith_mapped_grid))
    % Plot the base point
    view(3)
    % plot3(base_point(1,1),base_point(1,2),base_point(1,3),'r.','MarkerSize',50);
    quiver3(oriGIN(1,1),oriGIN(1,2),oriGIN(1,3), unit_normal_vectors(1,1),unit_normal_vectors(1,2),unit_normal_vectors(1,3),0,'g','Linewidth',3);
    quiver3(oriGIN(1,1),oriGIN(1,2),oriGIN(1,3), unit_vertical(1,1),unit_vertical(1,2),unit_vertical(1,3),0,'r','Linewidth',3);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
%     % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end