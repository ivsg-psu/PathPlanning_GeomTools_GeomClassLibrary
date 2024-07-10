%% 

% This code should be run section by section. 
% This code uses GPS class and geom class  
% This code is used to generate the plots for PPT
% load the data - easiest way: double click the .mat file and load the data
% Data Location: https://pennstateoffice365.sharepoint.com/:u:/r/sites/IntelligentVehiclesandSystemsGroup-Active/Shared%20Documents/IVSG/GitHubMirror/MappingVanDataCollection/ParsedData/2024-06-28/TestTrack_Entire_Loop.mat?csf=1&web=1&e=zIPfTp


% Revision History
% 
% Aneesh Batchu - 2024_07_09 
% -- wrote the code originally

%% Slide 1: Intro slide

% No plots here

%% Slide 2: Overview

% No plots here

%% Slide 3: Plot the Lidar Data in LLA coordinates 

% figure number
fig_num_LLA = 3000;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});
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
colors = z; 

figure(fig_num_slide_4);clf;
% Create the scatter plot
scatter3(x, y, z, 20, colors, 'filled');
colormap('jet'); 
colorbar; 
view(2); 
title('Z-coordinate is represented in colors');
xlabel('X[m]');
ylabel('Y[m]');

%% Slide 5: plot z-mean and angle difference in colors

fig_num_z_mean_slide_5 = 5000;
fig_num_angle_diff_slide_5 = 5001;

% LiDAR data
x = LiDAR_outer_edge(:,1);
y = LiDAR_outer_edge(:,2);
z = LiDAR_outer_edge(:,3); 

% Calculate the deviation of z-values from the mean
mean_z = mean(z);
deviation_from_mean = z - mean_z;

figure(fig_num_z_mean_slide_5); clf; 
% Create the scatter plot
scatter3(x, y, z, 20, deviation_from_mean, 'filled');
colormap('jet'); 
colorbar; 
view(2); 
title('Deviation of z-value from the mean');
xlabel('X[m]');
ylabel('Y[m]');

% plot angle between the unit vector of the data points and the vertical

% unit vectors of all the points
unit_vectors = fcn_geometry_calcUnitVector(LiDAR_outer_edge, -1);

% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_vectors.*unit_vector_vertical_direction,2);

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_unit_vectors_and_vertical = acos(dot_product);


% Calculate the angle difference from the vertical
vertical_angle = 90*(pi/180);
angle_difference = angle_btw_unit_unit_vectors_and_vertical - vertical_angle;

figure(fig_num_angle_diff_slide_5); clf;
% Create the scatter plot
scatter3(x, y, z, 20, angle_difference, 'filled');
colormap('jet'); 
colorbar; 
view(2); 
title('Angle between the unit vectors of the data points and the vertical');
xlabel('X[m]');
ylabel('Y[m]');

%% Slide 6: Plot hand-labeled boundary points on LLA

fig_num_LLA_boundary_points = 6000; 

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

% plot the points on LLA points
figure(fig_num_LLA_boundary_points)

xlabel('Latitude')
ylabel('Longitude')

% plot the LLA data
geoplot(LLA_data(:,1),LLA_data(:,2),'mo','MarkerSize',10);
hold on
geoplot(LLA_data(:,1),LLA_data(:,2),'k.','MarkerSize',10);

% plot the points in yellow color
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
% plot the same point in black over yellow points
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);

title('LLA Trace geometry with hand-labeled boundary points')

geobasemap satellite
geotickformat -dd

fig_bd_pts = 6001;
figure(fig_bd_pts);
clf; 

xlabel('Latitude')
ylabel('Longitude')

geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
hold on
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);

title('Boundary Points in LLA ')

geobasemap satellite
geotickformat -dd

%% Slide 7: Plot the boundary point in ENU coordinates on ENU LiDAR data

fig_num_slide_7 = 7000; 

% LiDAR data
x = LiDAR_outer_edge(:,1);
y = LiDAR_outer_edge(:,2);
z = LiDAR_outer_edge(:,3);

% Manually selected boundary points
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
ENU_data = gps_object.WGSLLA2ENU(boundary_points_LLA(:,1), boundary_points_LLA(:,2), boundary_points_LLA(:,3));

% Create the scatter plot
figure(fig_num_slide_7);clf;

% Height of the points
colors = z; 

scatter3(x, y, z, 20, colors, 'filled');
hold on
plot(ENU_data(:,1), ENU_data(:,2), 'k.', 'MarkerSize', 20);
colormap('jet'); 
colorbar; 
view(2); 

title('Boundary Points on ENU LiDAR data');
xlabel('X[m]');
ylabel('Y[m]');

%% Slide 8 and 9: Finding optimal value for point density based on mapped and unmapped grids

% This function generates the mapped and unmapped grids, and also drivable
% and non-drivable regions. 
% If fig_num = 1000, mapped and unmapped grids are plotted in figure(999), 
% and drivable and non-drivable in figure(1000)

% Figure number
fig_num = 8001;

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

figure(fig_num-1)
clf;
% plot the LiDAR data 
plot(input_points(:,1), input_points(:,2), 'b.', 'MarkerSize', 10);
view(2)
temp = axis;
%     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
axis_range_x = temp(2)-temp(1);
axis_range_y = temp(4)-temp(3);
percent_larger = 0.3;
axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% ------------------------------------------------------------------------
% Points per each grid: THE OPTIMAL VALUE OF THIS PARAMETER IS FOUND HERE
% ------------------------------------------------------------------------
point_density = 60;

% The maximum standard deviation limit of z_error (z - z_fit) after fitting a plane
std_threshold = 0.05;

% Maximum limit of the angle between vertical and the unit normal vector of
% the fitted plane
theta_threshold = 30*(pi/180); 

% Maximum and Minimum height threshold of mapped grids. The mean of Z
% values of mapped grids are calulated and compared with z_height_threshold
z_height_threshold = [mean(LiDAR_outer_edge(:,3))-0.5, mean(LiDAR_outer_edge(:,3))+0.5];
 
% flag for plotting in 3D
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z,...
    unmapped_gridlines,mapped_gridlines,drivable_gridlines,non_drivable_gridlines,gridCenters_unmapped_grids,mean_z_height_of_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, z_height_threshold, (flag_plot_in_3D), (fig_num));

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(angle_btw_unit_normals_and_vertical)>=1)
assert(length(standard_deviation_in_z)>=1)
assert(length(mean_z_height_of_mapped_grids)>=1)
assert(length(unmapped_gridlines)>=1)
assert(length(mapped_gridlines)>=1)
assert(length(drivable_gridlines)>=1)
assert(length(non_drivable_gridlines)>=1)
assert(length(gridCenters_unmapped_grids)>=1)

%% Slide 10: Find the optimal value of std_threshold for classifying mapped grids into drivable and non-drivable regions 

% In this slide, optimal value of std_threshold was determined

% This function generates the mapped and unmapped grids, and also drivable
% and non-drivable regions. 
% If fig_num = 1000, mapped and unmapped grids are plotted in figure(999), 
% and drivable and non-drivable in figure(1000)

% Figure number
fig_num = 10001;

% Flag for plotting hand-labeled boundary points
flag_plot_hand_labeled_boundary_points = 1; 

if flag_plot_hand_labeled_boundary_points
    % plot the ENU boundary points for reference
    figure(fig_num);
    clf;
    plot(ENU_data(:,1), ENU_data(:,2), 'k.', 'MarkerSize', 20);
end
% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.265; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% Points per each grid
point_density = 60;

% ------------------------------------------------------------------------
% The maximum standard deviation limit of z_error (z - z_fit) after fitting a plane: THE OPTIMAL VALUE OF THIS PARAMETER IS FOUND HERE
% ------------------------------------------------------------------------
std_threshold = 0.03; 

% Maximum limit of the angle between vertical and the unit normal vector of
% the fitted plane
theta_threshold =[]; %30*(pi/180); 

% Maximum and Minimum height threshold of mapped grids. The mean of Z
% values of mapped grids are calulated and compared with z_height_threshold
z_height_threshold = []; %[mean(LiDAR_outer_edge(:,3))-0.5, mean(LiDAR_outer_edge(:,3))+0.5];
 
% flag for plotting in 3D
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z,...
    unmapped_gridlines,mapped_gridlines,drivable_gridlines,non_drivable_gridlines,gridCenters_unmapped_grids,mean_z_height_of_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, z_height_threshold, (flag_plot_in_3D), (fig_num));

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(isempty(angle_btw_unit_normals_and_vertical))
assert(length(standard_deviation_in_z)>=1)
assert(length(mean_z_height_of_mapped_grids)>=1)
assert(length(unmapped_gridlines)>=1)
assert(length(mapped_gridlines)>=1)
assert(length(drivable_gridlines)>=1)
assert(length(non_drivable_gridlines)>=1)
assert(length(gridCenters_unmapped_grids)>=1)

%% Slide 11: Find the optimal value of theta_threshold for classifying mapped grids into drivable and non-drivable regions 

% In this slide, optimal value of theta_threshold was determined

% This function generates the mapped and unmapped grids, and also drivable
% and non-drivable regions. 
% If fig_num = 1000, mapped and unmapped grids are plotted in figure(999), 
% and drivable and non-drivable in figure(1000)

% Figure number
fig_num = 11001;

% Flag for plotting hand-labeled boundary points
flag_plot_hand_labeled_boundary_points = 1; 

if flag_plot_hand_labeled_boundary_points
    % plot the ENU boundary points for reference
    figure(fig_num);
    clf;
    plot(ENU_data(:,1), ENU_data(:,2), 'k.', 'MarkerSize', 20);
end

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.265; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% Points per each grid
point_density = 60;

% The maximum standard deviation limit of z_error (z - z_fit) after fitting a plane
std_threshold = []; %0.05; 

% ------------------------------------------------------------------------
% Maximum limit of the angle between vertical and the unit normal vector of
% the fitted plane:                                                        THE OPTIMAL VALUE OF THIS PARAMETER IS FOUND HERE
% ------------------------------------------------------------------------
theta_threshold =50*(pi/180); 

% Maximum and Minimum height threshold of mapped grids. The mean of Z
% values of mapped grids are calulated and compared with z_height_threshold
z_height_threshold = [];%[mean(LiDAR_outer_edge(:,3))-0.5, mean(LiDAR_outer_edge(:,3))+0.5];
 
% flag for plotting in 3D
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z,...
    unmapped_gridlines,mapped_gridlines,drivable_gridlines,non_drivable_gridlines,gridCenters_unmapped_grids,mean_z_height_of_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, z_height_threshold, (flag_plot_in_3D), (fig_num));

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(angle_btw_unit_normals_and_vertical)>=1)
assert(length(standard_deviation_in_z)>=1)
assert(length(mean_z_height_of_mapped_grids)>=1)
assert(length(unmapped_gridlines)>=1)
assert(length(mapped_gridlines)>=1)
assert(length(drivable_gridlines)>=1)
assert(length(non_drivable_gridlines)>=1)
assert(length(gridCenters_unmapped_grids)>=1)

%% Slide 12: Find the optimal value of z_height_threshold for classifying mapped grids into drivable and non-drivable regions 

% In this slide, optimal value of z_height_threshold was determined

% This function generates the mapped and unmapped grids, and also drivable
% and non-drivable regions. 
% If fig_num = 1000, mapped and unmapped grids are plotted in figure(999), 
% and drivable and non-drivable in figure(1000)

% Figure number
fig_num = 12001;

% Flag for plotting hand-labeled boundary points
flag_plot_hand_labeled_boundary_points = 1; 

if flag_plot_hand_labeled_boundary_points
    % plot the ENU boundary points for reference
    figure(fig_num);
    clf;
    plot(ENU_data(:,1), ENU_data(:,2), 'k.', 'MarkerSize', 20);
end

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.265; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% Points per each grid
point_density = 60;

% The maximum standard deviation limit of z_error (z - z_fit) after fitting a plane
std_threshold = []; %0.05; 


% Maximum limit of the angle between vertical and the unit normal vector of
% the fitted plane:                                                        
theta_threshold =[]; %70*(pi/180); 

% ------------------------------------------------------------------------
% Maximum and Minimum height threshold of mapped grids. The mean of Z
% values of mapped grids are calulated and compared with z_height_threshold:THE OPTIMAL VALUE OF THIS PARAMETER IS FOUND HERE
% ------------------------------------------------------------------------
z_height_threshold = [mean(LiDAR_outer_edge(:,3))-0.5, mean(LiDAR_outer_edge(:,3))+0.5];
 
% flag for plotting in 3D
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z,...
    unmapped_gridlines,mapped_gridlines,drivable_gridlines,non_drivable_gridlines,gridCenters_unmapped_grids,mean_z_height_of_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, z_height_threshold, (flag_plot_in_3D), (fig_num));

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
assert(length(unmapped_grids)>=1)
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(isempty(angle_btw_unit_normals_and_vertical))
assert(length(standard_deviation_in_z)>=1)
assert(length(mean_z_height_of_mapped_grids)>=1)
assert(length(unmapped_gridlines)>=1)
assert(length(mapped_gridlines)>=1)
assert(length(drivable_gridlines)>=1)
assert(length(non_drivable_gridlines)>=1)
assert(length(gridCenters_unmapped_grids)>=1)

%% Slide 13 and 14: Find the boundary points of drivable and non-drivable grids by choosing the optimal values for thresholds. 

% This function generates the mapped and unmapped grids, and also drivable
% and non-drivable regions. 
% If fig_num = 1000, mapped and unmapped grids are plotted in figure(999), 
% and drivable and non-drivable in figure(1000)

% Figure number for plotting grids
fig_num = 12001;

% Figure numbers for boundary points plots
fig_num_drivable_non_drivable = 700;
fig_num_mapped_unmapped = 701; 
fig_num_bd_pts_ENU = 702; 

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.265; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% Points per each grid
point_density = 60;

% The maximum standard deviation limit of z_error (z - z_fit) after fitting a plane
std_threshold = 0.05; 

% Maximum limit of the angle between vertical and the unit normal vector of
% the fitted plane:                                                        
theta_threshold = [];%70*(pi/180); 

% Maximum and Minimum height threshold of mapped grids. The mean of Z
% values of mapped grids are calulated and compared with z_height_threshold
z_height_threshold = [];%[mean(LiDAR_outer_edge(:,3))-0.5, mean(LiDAR_outer_edge(:,3))+0.5];
 
% flag for plotting in 3D
flag_plot_in_3D = 0; 

[drivable_grids,non_drivable_grids,unmapped_grids,gridCenters_mapped_grids,drivable_grid_numbers_in_mapped_grids,...
    non_drivable_grid_numbers_in_mapped_grids,angle_btw_unit_normals_and_vertical,standard_deviation_in_z,...
    unmapped_gridlines,mapped_gridlines,drivable_gridlines,non_drivable_gridlines,gridCenters_unmapped_grids,mean_z_height_of_mapped_grids] = fcn_geometry_surfaceAnalysis(input_points, grid_size, grid_boundaries, point_density, theta_threshold, std_threshold, z_height_threshold, (flag_plot_in_3D), (fig_num));

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)>=1)
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(angle_btw_unit_normals_and_vertical)>=1)
% assert(length(standard_deviation_in_z)>=1)
% assert(length(mean_z_height_of_mapped_grids)>=1)
% assert(length(unmapped_gridlines)>=1)
% assert(length(mapped_gridlines)>=1)
% assert(length(drivable_gridlines)>=1)
% assert(length(non_drivable_gridlines)>=1)
% assert(length(gridCenters_unmapped_grids)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of mapped and unmapped grids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XYZ_matrix_mapped_grids = [gridCenters_mapped_grids(:,1:2) ones(length(gridCenters_mapped_grids),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_unmapped_grids(:,1:2) gridCenters_unmapped_grids(:,4)]; 

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

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_mapped_unmapped);
clf;
boundary_points_mapped_unmapped = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)==zeros(0,1))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)]; 

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
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);

assert(length(drivable_grids)>=1)
assert(length(non_drivable_grids)>=1)
% assert(length(unmapped_grids)==zeros(0,1))
assert(length(gridCenters_mapped_grids(:,1))>=1)
assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[members, id_x] = ismember(boundary_points,boundary_points_mapped_unmapped,'rows'); 

not_boundary_points = boundary_points(members,:);

true_boundary_points = boundary_points(members==0,:);

% figure(fig_num_bd_pts_ENU)
% clf;
figure(10001)
hold on
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40)
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30)


 %% Slide 13 and 14: Plot computed boundary points and hand-labled in LLA

fig_final_bd_pts = 113; 

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


Trace_coordinates = [true_boundary_points,zeros(length(true_boundary_points),1)]; 

 % Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_data_computed_boundary_pts = gps_object.ENU2WGSLLA(Trace_coordinates);

% Plot the LLA boubndary points
figure(fig_final_bd_pts);
clf

xlabel('Latitude')
ylabel('Longitude')

geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'y.','MarkerSize',30);
hold on
geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);

geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'c.','MarkerSize',40);
geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'b.','MarkerSize',30);

title('Boundary Points in LLA ')

geobasemap satellite
geotickformat -dd