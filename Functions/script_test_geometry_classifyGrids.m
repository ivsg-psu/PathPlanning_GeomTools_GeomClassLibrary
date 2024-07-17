
%% Get the LiDAR data based on ring

% INPUTS: LiDAR_Scan_ENU_Entire_Loop, ringID, row_start, row_end, fig_num (use varargin)
% OUTPUTS: LiDAR_allPoints_of_ring
% Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;


% figure number
fig_num_LLA = 3001;
ring = 15; 

% Fill raw LiDAR data
row_start = 1400;
row_end   = 1410;
LiDAR_scan_line = [];
concatenate_LiDAR_points = [];
LiDAR_ring = []; 
for ith_row = row_start:row_end
    LiDAR_ring = [LiDAR_ring; LiDAR_ENU_cell{ith_row}(:,5); nan*LiDAR_ENU_cell{ith_row}(1,5)]; %#ok<AGROW>
    LiDAR_scan_line   = [LiDAR_scan_line; (ith_row)*ones(length(LiDAR_ENU_cell{ith_row}(:,1:3))+1,1)];%#ok<AGROW>
   concatenate_LiDAR_points = [concatenate_LiDAR_points; LiDAR_ENU_cell{ith_row}(:,1:3); nan*LiDAR_ENU_cell{ith_row}(1,1:3)]; %#ok<AGROW>
end

% LiDAR data: contains scan line number, ENU data, LiDAR ring
LiDAR_allPoints = [LiDAR_scan_line,concatenate_LiDAR_points,LiDAR_ring];

% Select the data based on the ring

LiDAR_index_ring_indices = find(LiDAR_allPoints(:,5) == ring);  

% LiDAR points of the ring
LiDAR_allPoints_of_ring = LiDAR_allPoints(LiDAR_index_ring_indices,:); 


% % Define GPS object
% reference_latitude = 40.86368573;
% reference_longitude = -77.83592832;
% reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_LIDAR_allPoints_of_ring = gps_object.ENU2WGSLLA(LiDAR_allPoints_of_ring(:,2:4));


% Plot the ENU results
figure(fig_num_LLA);clf;

xlabel('Latitude')
ylabel('Longitude')

geoplot(LLA_LIDAR_allPoints_of_ring(:,1),LLA_LIDAR_allPoints_of_ring(:,2),'mo','MarkerSize',10);
hold on
geoplot(LLA_LIDAR_allPoints_of_ring(:,1),LLA_LIDAR_allPoints_of_ring(:,2),'k.','MarkerSize',10);

title('LLA Trace geometry')

geobasemap satellite

%% Gather all the raw and drivable data and pot



%% INPUTS and Plots the LLA data

%% Find the drivable area, approximately. Each cell array is ENU XYZ data.

LiDAR_allPoints_ENU = LiDAR_allPoints_of_ring(:,2:4); 

test_scan = LiDAR_allPoints_ENU;

figure(233);
plot3(test_scan(:,1), test_scan(:,2), test_scan(:,3), 'k.', 'MarkerSize',20)
figure(34545);
clf;
% The following numbers were found manually and depend on the LIDAR being
% used
drivable_row_start = 400;
drivable_row_end   = 743;

subplot(1,2,1)
hold on;
grid on;
axis equal;

plot(test_scan(:,1),test_scan(:,3),'k.','MarkerSize',20);
plot(test_scan(drivable_row_start:drivable_row_end,1),...
    test_scan(drivable_row_start:drivable_row_end,3),'g.','MarkerSize',10);


subplot(1,2,2)
hold on;
grid on;
axis equal;

plot(test_scan(:,2),test_scan(:,3),'k.','MarkerSize',20);
plot(test_scan(drivable_row_start:drivable_row_end,2),...
    test_scan(drivable_row_start:drivable_row_end,3),'g.','MarkerSize',10);


%%


% Fill drivable area LiDAR data 
LiDAR_drivablePoints = [];
for ith_row = row_start:row_end
    LiDAR_drivablePoints = [LiDAR_drivablePoints; ...
        LiDAR_ENU_cell{ith_row}(drivable_row_start:drivable_row_end,:); nan*LiDAR_ENU_cell{ith_row}(1,:)]; %#ok<AGROW>
end


% LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});

% Input points (LiDAR data)
input_points = LiDAR_allPoints(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

grid_boundaries = [Min_x Max_x Min_y Max_y]; 

% As of now. 
% point_density = 20;

% LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{:});
Trace_coordinates = LiDAR_allPoints;

% Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_LIDAR_allPoints = gps_object.ENU2WGSLLA(Trace_coordinates);
LLA_LIDAR_drivablePoints = gps_object.ENU2WGSLLA(LiDAR_drivablePoints);


% Plot the ENU results
figure(fig_num_LLA);clf;

xlabel('Latitude')
ylabel('Longitude')

geoplot(LLA_LIDAR_allPoints(:,1),LLA_LIDAR_allPoints(:,2),'mo','MarkerSize',10);
hold on
geoplot(LLA_LIDAR_allPoints(:,1),LLA_LIDAR_allPoints(:,2),'k.','MarkerSize',10);
geoplot(LLA_LIDAR_drivablePoints(:,1),LLA_LIDAR_drivablePoints(:,2),'g.','MarkerSize',10);

title('LLA Trace geometry')

geobasemap satellite

%% Hand Labeled Boundary points

boundary_points_LLA = [40.865759464410559, -77.830964010633224, 0;
    40.865755710442414, -77.830956939857998, 0;
    40.865752421719513, -77.830948769793835, 0;
    40.865748357717237, -77.830941908178062, 0;
    40.865745581471145, -77.830933910752478, 0;
    40.865740732632631, -77.830927575313964, 0;
    40.865736995531826, -77.830919982286616, 0;
    40.865734153848372, -77.830912003671742, 0;
    40.865730068320772, -77.830905010048482, 0;
    40.865725645797589, -77.830897863602658, 0;
    40.865722297500483, -77.830890350824461, 0];

% figure number
fig_num_LLA = 3001;
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
hold on
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);
%%

fig_bd_pts = 3300;
figure(fig_bd_pts);

xlabel('Latitude')
ylabel('Longitude')

geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
hold on
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);

title('Boundary Points in LLA ')

geobasemap satellite
geotickformat -dd
%% - Plot drivable path of 1400:1409


drivable_path_1400_1410 =  [40.865704401971719 -77.830904080921826 0;
    40.865708684754139 -77.830910955129141 0;
    40.865712687063684 -77.830917925434548 0;
    40.865716658256282 -77.830925002208645 0;
    40.865720707785179 -77.830932010671759 0;
    40.865724663481600 -77.830939007034146 0;
    40.865728496090632 -77.830946206951452 0;
    40.865732318250444 -77.830953304121806 0;
    40.865736101077680 -77.830960515870203 0;
    40.865739763290073 -77.830967819269759 0;
    40.865743204622873 -77.830975304540672 0;
    40.865682626083093 -77.830927132797171 0;
    40.865686683570189 -77.830933708132505 0;
    40.865690680995193 -77.830940265372661 0;
    40.865694249441972 -77.830947400549462 0;
    40.865698048099269 -77.830954161967071 0;
    40.865702011711655 -77.830960627038735 0;
    40.865705536069854 -77.830967669676696 0;
    40.865709324637606 -77.830974538222733 0;
    40.865713105884005 -77.830981387897069 0;
    40.865716549822061 -77.830988593573522 0;
    40.865720548721377 -77.830995199658403 0];

figure(fig_num_LLA)
geoplot(drivable_path_1400_1410(:,1),drivable_path_1400_1410(:,2),'g.','MarkerSize',40);
hold on
geoplot(drivable_path_1400_1410(:,1),drivable_path_1400_1410(:,2),'k.','MarkerSize',30);

%%

fig_num_slide_7 = 7000; 

% LiDAR data
x = LiDAR_allPoints(:,1);
y = LiDAR_allPoints(:,2);
z = LiDAR_allPoints(:,3);

gps_object = GPS(); % Initiate the class object for GPS

% Use the class to convert LLA to ENU
ENU_drivable_path_1400_1410 = gps_object.WGSLLA2ENU(drivable_path_1400_1410(:,1), drivable_path_1400_1410(:,2), drivable_path_1400_1410(:,3));

% Create the scatter plot
figure(fig_num_slide_7);clf;

% Height of the points
colors = z; 

scatter3(x, y, z, 20, colors, 'filled');
hold on
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);
colormap('jet'); 
colorbar; 
view(2); 

title('Boundary Points on ENU LiDAR data');
xlabel('X[m]');
ylabel('Y[m]');

%% MAIN - grids with zero and non-zero

% INPUTS - input_points,grid_size,grid_boundaries, plotting
% OUTPUTS - gridIndices_cell_array, 
% total_N_points_in_each_grid,gridCenters, grids_with_zero_points, grids_greater_than_zero_points 
% fcn_geometry_findGridsWithPoints


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
gridCenters_zero_point_density = gridCenters(grids_with_zero_points,1:2); 

% Find grids with more than zero points
grids_greater_than_zero_points = find((total_N_points_in_each_grid(:,1) > 0)); 

% Grid Centers of the grids with zero point density (Unmapped grid centers)
gridCenters_greater_than_zero_point_density = gridCenters(grids_greater_than_zero_points,1:2); 

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
plot_gridCenters_zero_density = [gridCenters_zero_point_density, zeros(length(gridCenters_zero_point_density(:,1)),1)];
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_zero_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids greater than zero point density';
legend_position = [];
plot_gridCenters_greater_than_zero_density = [gridCenters_greater_than_zero_point_density, zeros(length(gridCenters_greater_than_zero_point_density(:,1)),1)];
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_greater_than_zero_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% Plot mapped grid lines

figure(7890);clf;
hold on 
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the mapped grids')
length_points_in_domain_grids_greater_than_zero = zeros(length(grids_greater_than_zero_points),1); 
mapped_gridlines = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11
    for ith_domain = 1:length(grids_greater_than_zero_points)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

        current_color = [0.2 0.2 0.2];
        % Plot current AABB
        current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_domain),1:4);

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
        rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_domain);
        points_in_domain = input_points(rows_in_domain,:);

        length_points_in_domain_grids_greater_than_zero(ith_domain) = length(points_in_domain);

        length_gridlines = length(gridlines);
        % % Plot the mapped points green
        % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

        % Plot the result
        % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        % hold on
        plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
        mapped_gridlines(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;
    end


    for ith_text = 1:length(grids_greater_than_zero_points(:,1))
        current_text = sprintf('%.0d',ith_text);
        % Place the text on the grid center
        text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
    end
    % plot the drivable path
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

    point_density = median(length_points_in_domain_grids_greater_than_zero)/2; 

    % point_density = 20; 
%% MAIN - classify grids with more than zero points into mapped and unmapped


% INPUTS -
% point_density,total_N_points_in_each_grid,grids_greater_than_zero_points,
% grid_centers, plotting options
% OUTPUTS - original_grids_with_required_point_density
% fcn_geometry_classifyGridsintoMappedUnmapped

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
gridCenters_low_point_density = gridCenters(original_grids_with_low_point_density,1:2);  

% Grid Centers of the grids with required point density (Mapped grid centers)
gridCenters_required_point_density = gridCenters(original_grids_with_required_point_density,1:2);  

% Current grid numbers of the grids with low point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_low_point_density = find((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Current grid numbers of the grids with required point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_required_point_density = find(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% plotting
fig_num = 121; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Unmapped grids';
legend_position = [];
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_low_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Mapped grids';
legend_position = [];
plot_gridCenters_required_point_density = [gridCenters_required_point_density, zeros(length(gridCenters_required_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_required_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% ENU - plotting mapped and unmapped
figure(1234)
clf;
plot(gridCenters_low_point_density(:,1), gridCenters_low_point_density(:,2), '.','MarkerSize',40,'Color',[0.8 0.8 0.8]);
hold on
grid on
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

% plot the drivable path
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

%% MAIN - classify mapped grids into drivable and non-drivable

% INPUT
% original_grids_with_required_point_density
input_points = LiDAR_allPoints; 
std_threshold = []; 
theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters


% OUTPUTS
% standard_deviation_in_z,angle_btw_unit_normals_and_vertical,original_drivable_grids,original_non_drivable_grids  
% current_drivable_grid_numbers_in_mapped_grids,current_non_drivable_grid_numbers_in_mapped_grids,gridCenters_drivable_grids,gridCenters_non_drivable_grids

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_grids_with_required_point_density); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


% z_height of all the points 
% mean_z_of_mapped_grids = zeros(total_mapped_grids,1); 

% Loop through all the mapped grids, recording standard deviation, unit
% vectors 
% if 0==flag_max_speed
%     h_waitbar = waitbar(0,'Performing surface analysis...');
% end

for ith_mapped_grid = 1:total_mapped_grids
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% standard_deviation_in_z = round(standard_deviation_in_z,4); 
if ~isempty(std_threshold) && isempty(theta_threshold)

    % STEP 1: Standard deviation of the orthogonal (perpendicular) distances of
    % the points to the plane (after fit)
    % Find the grids that are within standard deviation limit
    % This is not enough (delta Y) is also important
    % mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < std_threshold;
    mapped_grids_within_std_threshold = standard_deviation_in_z < std_threshold;
    
    % Grids that satisy the conditions of (STEP 1). The grids that
    % are within the std threshold
    mapped_grids_within_all_thresholds = (mapped_grids_within_std_threshold == 1);

    % The angle between unit vertical and the unit_normal_vector is computed to
    % determine how close the normal vector is to vertical direction. In
    % this case, the angle between unit normals and vertical is empty. 
    angle_btw_unit_normals_and_vertical = [];

elseif isempty(std_threshold) && ~isempty(theta_threshold)

    % STEP 2
    % Comparing normal vector with vertical direction
    unit_vector_vertical_direction = [0 0 1];

    % The dot product is computed to find the angle between the vectors
    dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);
  
    % The angle between unit vertical and the unit_normal_vector is computed to
    % determine how close the normal vector is to vertical direction.
    angle_btw_unit_normals_and_vertical = acos(dot_product);

    % Find the grids (with a fitted plane) that are within the vertical
    % threshold (change to a different name: vertical threshold)
    mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < theta_threshold;

    % Grids that satisy the conditions (STEP 2). The grids that
    % are within the vertical threshold
    mapped_grids_within_all_thresholds = (mapped_grids_within_vertical_threshold == 1);

else

    % STEP 1: Standard deviation of the orthogonal (perpendicular) distances of
    % the points to the plane (after fit)
    % Find the grids that are within standard deviation limit
    % This is not enough (delta Y) is also important
    % mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < std_threshold;
    mapped_grids_within_std_threshold = standard_deviation_in_z < std_threshold;

    % STEP 2
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

    % Grids that satisy the conditions of (STEP 1 & STEP 2). The grids that
    % are within the standar deviation and vertical threshold
    mapped_grids_within_vertical_and_std_thresholds = (mapped_grids_within_vertical_threshold == 1) & (mapped_grids_within_std_threshold == 1);
    
    % mapped grids within all the thresholds 
    mapped_grids_within_all_thresholds = mapped_grids_within_vertical_and_std_thresholds;

end

% Find the drivable grids (original)
original_drivable_grids = original_grids_with_required_point_density(mapped_grids_within_all_thresholds); 

% Find the non-drivable grids (original)
original_non_drivable_grids = original_grids_with_required_point_density(mapped_grids_within_all_thresholds == 0);

% Final drivable grid numbers of the mapped grids
current_drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, original_drivable_grids));

% Final non drivable grid numbers of the mapped grids
current_non_drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, original_non_drivable_grids));

% Grid centers of drivable grids 
gridCenters_drivable_grids = [gridCenters(original_drivable_grids,1), gridCenters(original_drivable_grids,2), ones(length(original_drivable_grids),1)]; 

% Grid centers of nondrivable grids
gridCenters_non_drivable_grids = [gridCenters(original_non_drivable_grids,1), gridCenters(original_non_drivable_grids,2), zeros(length(original_non_drivable_grids),1)]; 

% Concatenate the grid centers of drivable and non-drivable grids (2D)
gridCenters_mapped_grids = [gridCenters_drivable_grids; gridCenters_non_drivable_grids];

% plotting
fig_num = 6000; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0, 1, 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [1, 0, 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
plot_gridCenters_non_drivable_grids = [gridCenters_non_drivable_grids(:,1:2), zeros(length(gridCenters_non_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

%% Boundary points on LLA

boundary_points_LLA = [40.865759464410559, -77.830964010633224, 0;
    40.865755710442414, -77.830956939857998, 0;
    40.865752421719513, -77.830948769793835, 0;
    40.865748357717237, -77.830941908178062, 0;
    40.865745581471145, -77.830933910752478, 0;
    40.865740732632631, -77.830927575313964, 0;
    40.865736995531826, -77.830919982286616, 0;
    40.865734153848372, -77.830912003671742, 0;
    40.865730068320772, -77.830905010048482, 0;
    40.865725645797589, -77.830897863602658, 0;
    40.865722297500483, -77.830890350824461, 0];

gps_object = GPS(); % Initiate the class object for GPS

% Use the class to convert LLA to ENU
ENU_boundary_points_data = gps_object.WGSLLA2ENU(boundary_points_LLA(:,1), boundary_points_LLA(:,2), boundary_points_LLA(:,3));


% plot handlabeled boundar points
marker_size = 30;
RGB_triplet = [1, 1, 0]; 
legend_option = 1;
legend_name = 'Hand-labeled boundary points';
legend_position = [];

[~] = fcn_geometry_plotGridCenters(ENU_boundary_points_data,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);


% geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
% geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);


%% Boundary Points

%% Pre the grid centers for finding boundary points

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_mapped_unmapped = 767787; 

XYZ_matrix_mapped_grids = [gridCenters_required_point_density(:,1:2) ones(length(gridCenters_required_point_density(:,1)),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_low_point_density(:,1:2) zeros(length(gridCenters_low_point_density(:,1)),1)]; 

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

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)==zeros(0,1))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_num_drivable_non_drivable = 98898;

% XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)];
XYZ_matrix = gridCenters_mapped_grids;

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

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)==zeros(0,1))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% % True boundary points
% Trace_coordinates = [true_boundary_points,zeros(length(true_boundary_points),1)]; 
% 
%  % Define GPS object
% reference_latitude = 40.86368573;
% reference_longitude = -77.83592832;
% reference_altitude = 344.189;
% gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class
% 
% % Use the class to convert LLA to ENU
% LLA_data_computed_boundary_pts = gps_object.ENU2WGSLLA(Trace_coordinates);
% 
% 
% figure(6000)
% 
% geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'c.','MarkerSize',30);
% geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'b.','MarkerSize',15);


%% ENU mapped grids

figure(12345)
clf;
hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of drivable and non-drivable grids in ENU')
% plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);


mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1)); 
for ith_text = 1:length(gridCenters_required_point_density(:,1))
    current_text = sprintf('%.0d',mapped_grid_numbers(ith_text));
     % Place the text on the grid center
    text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold');
end

% % plot the drivable path
% plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 30);
% plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 10);

%% ENU mapped grids - with std_deviation of z in grids as the text 

% figure(12345)
% clf;
% hold on
% axis on
% grid on 
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Grid centers of the mapped grids')
% plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% standard_deviation_in_z_round = round(standard_deviation_in_z,3);
% % mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1)); 
% for ith_text = 1:length(gridCenters_required_point_density(:,1))
%     current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
%      % Place the text on the grid center
%     text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[0.4660 0.6740 0.1880],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
% end
% 
% % plot the drivable path
% plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
% plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);
% 
% temp = axis;
% %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
% axis_range_x = temp(2)-temp(1);
% axis_range_y = temp(4)-temp(3);
% percent_larger = 0.5;
% axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
% 


fig_num = 34332; 
figure(fig_num)

hold on 
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')
length_points_in_domain_grids_required_density = zeros(length(original_grids_with_required_point_density),1); 
mapped_gridlines = zeros(11*length(original_grids_with_required_point_density),2); % length(gridlines) = 11
    for ith_domain = 1:length(original_grids_with_required_point_density)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

        current_color = [0.2 0.2 0.2];
        % Plot current AABB
        current_AABB = grid_AABBs(original_grids_with_required_point_density(ith_domain),1:4);

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
        rows_in_domain = gridIndices==original_grids_with_required_point_density(ith_domain);
        points_in_domain = input_points(rows_in_domain,:);

        length_points_in_domain_grids_required_density(ith_domain) = length(points_in_domain);

        length_gridlines = length(gridlines);
        % % Plot the mapped points green
        % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

        % Plot the result
        % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        % hold on
        % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
        mapped_gridlines(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;
    end

    standard_deviation_in_z_round = round(standard_deviation_in_z,3);
    % mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
    for ith_text = 1:length(gridCenters_required_point_density(:,1))
        current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
        % Place the text on the grid center
        text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[0.4660 0.6740 0.1880],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
    end


    % plot the drivable path
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

%% ENU mapped grids - with angle_btw_unit_normals_and_vertical in all grids as the text 

figure(1277)
clf;
hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Angle between unit normal and vertical on Grid centers of the mapped grids')
plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
angle_btw_unit_normals_and_vertical_round = round((angle_btw_unit_normals_and_vertical*180/pi),3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1)); 
for ith_text = 1:length(gridCenters_required_point_density(:,1))
    current_text = sprintf('%.0f',angle_btw_unit_normals_and_vertical_round(ith_text));
     % Place the text on the grid center
    text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[0.4660 0.6740 0.1880],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end

% plot the drivable path
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

%%

fig_num = 34337; 
figure(fig_num)

hold on 
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')
length_points_in_domain_grids_required_density = zeros(length(original_grids_with_required_point_density),1); 
mapped_gridlines = zeros(11*length(original_grids_with_required_point_density),2); % length(gridlines) = 11
    for ith_domain = 1:length(original_grids_with_required_point_density)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

        current_color = [0.2 0.2 0.2];
        % Plot current AABB
        current_AABB = grid_AABBs(original_grids_with_required_point_density(ith_domain),1:4);

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
        rows_in_domain = gridIndices==original_grids_with_required_point_density(ith_domain);
        points_in_domain = input_points(rows_in_domain,:);

        length_points_in_domain_grids_required_density(ith_domain) = length(points_in_domain);

        length_gridlines = length(gridlines);
        % % Plot the mapped points green
        % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

        % Plot the result
        % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        % hold on
        % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
        mapped_gridlines(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;
    end

    angle_btw_unit_normals_and_vertical_round = round((angle_btw_unit_normals_and_vertical*180/pi),3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1)); 
for ith_text = 1:length(gridCenters_required_point_density(:,1))
    current_text = sprintf('%.0f',angle_btw_unit_normals_and_vertical_round(ith_text));
     % Place the text on the grid center
    text(gridCenters_required_point_density(ith_text,1), gridCenters_required_point_density(ith_text,2),current_text,'Color',[0.4660 0.6740 0.1880],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end

% plot the drivable path
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

%%


fig_num = 34367; 
figure(fig_num)


% input_points = LiDAR_outer_edge; 
fig_num_slide_4 = fig_num; 

% LiDAR data
x = LiDAR_allPoints(:,1);
y = LiDAR_allPoints(:,2);
z = LiDAR_allPoints(:,3);

% Height of the points
colors = z; 

figure(fig_num_slide_4);
% Create the scatter plot
scatter3(x, y, z, 20, colors, 'filled');
colormap('jet'); 
colorbar; 
view(2); 
title('Z-coordinate is represented in colors');
xlabel('X[m]');
ylabel('Y[m]');

hold on 
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')
length_points_in_domain_grids_required_density = zeros(length(original_non_drivable_grids),1); 
mapped_gridlines = zeros(11*length(original_non_drivable_grids),2); % length(gridlines) = 11
    for ith_domain = 1:length(original_non_drivable_grids)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

        current_color = [0.2 0.2 0.2];
        % Plot current AABB
        current_AABB = grid_AABBs(original_non_drivable_grids(ith_domain),1:4);

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
        rows_in_domain = gridIndices==original_non_drivable_grids(ith_domain);
        points_in_domain = input_points(rows_in_domain,:);

        length_points_in_domain_grids_required_density(ith_domain) = length(points_in_domain);

        length_gridlines = length(gridlines);
        % % Plot the mapped points green
        % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

        % Plot the result
        % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        % hold on
        % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
        mapped_gridlines(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;
    end

    % plot the drivable path
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'g.', 'MarkerSize', 40);
    plot(ENU_drivable_path_1400_1410(:,1), ENU_drivable_path_1400_1410(:,2), 'k.', 'MarkerSize', 20);

%%
oriGIN = [0 0 0]; 
unit_vertical = [0 0 1]; 
input_points = LiDAR_allPoints;
total_non_drivable_grids = length(original_non_drivable_grids);
for ith_mapped_grid = 75
    figure(ith_mapped_grid*10000); clf;

    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors, base_point, ~] =...
        fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),ith_mapped_grid*10000);
    xlabel('X[m]')
    ylabel('Y[m]')
    zlabel('Z[m]')
    title(sprintf('Plane fit of  %dth mapped grid',ith_mapped_grid))
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
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end




