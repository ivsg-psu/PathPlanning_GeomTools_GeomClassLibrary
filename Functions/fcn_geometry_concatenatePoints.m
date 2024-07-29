function [concatenate_LiDAR_XYZ_points,concatenate_VehiclePose_XYZ_points...
    ,concatenate_unit_ortho_vehicle_vectors_XYZ,intensity_min, intensity_max] = ...
    fcn_geometry_concatenatePoints(LiDAR_Scan_ENU_Entire_Loop,VehiclePose,...
    rings_to_analyze,scanLineNumber_start,scanLineNumber_end,unit_ortho_vehicle_vectors_XY,varargin)
%% fcn_geometry_concatenatePoints
% Takes in the inputs of where on the track we want to analyze and provides
% a concatenated list of points for that area. 
%
% INPUTS:
% LiDAR_Scan_ENU_Entire_Loop
% VehiclePose
% rings_to_analyze
% scanLineNumber_start
% scanLineNumber_end
% unit_ortho_vehicle_vectors_XY
% 
% OPTIONAL INPUTS
% (fig_num)
% (ENU_3D_fig_num)
% (simple_view_toggle_1)
% (scaling_1)
% (ENU_XZ_fig_num)
% (ENU_YZ_fig_num)
% (scaling_2)
% (reference_longitude)
% (reference_latitude)
% (reference_altitude)
% (LLA_fig_num)
%
% REQUIREMENTS: 
% fcn_geometry_vehicleOrientation for unit_ortho_vehicle_vectors_XY
%
% EXAMPLES:
% script_test_fcn_geometry_concatenatePoints
%
% 2024_07_17 - S Brennan
% -- wrote the code
% 2024_07_24 - A Goncharov
% -- Functionalized this part of the code
%

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==17 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(6,17);
 
    end
end

%ENU_3D_fig_num 
fig_num = [];
if (7<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
    end
end

%ENU_3D_fig_num 
ENU_3D_fig_num = [];
if (8<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        ENU_3D_fig_num = temp;
    end
end

%simple view toggle
simple_view_toggle_1 = 1;
if (9<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        simple_view_toggle_1 = temp;
    end
end

%Scaling function, default is 3
scaling_1 = 3;
if (10<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        scaling_1 = temp;
    end
end

ENU_XZ_fig_num = [];
if (11<=nargin)
    temp = varargin{5};
    if ~isempty(temp)
        ENU_XZ_fig_num = temp;
    end
end

ENU_YZ_fig_num = [];
if (12<=nargin)
    temp = varargin{6};
    if ~isempty(temp)
        ENU_YZ_fig_num = temp;
    end
end

scaling_2 = 3;
if (13<=nargin)
    temp = varargin{7};
    if ~isempty(temp)
        scaling_2 = temp;
    end
end


reference_latitude = 40.86368573;
if (14<=nargin)
    temp = varargin{8};
    if ~isempty(temp)
        reference_latitude = temp;
    end
end

reference_longitude = -77.83592832;
if (13<=nargin)
    temp = varargin{9};
    if ~isempty(temp)
        reference_longitude = temp;
    end
end

reference_altitude = 344.189;

if (13<=nargin)
    temp = varargin{10};
    if ~isempty(temp)
        reference_altitude = temp;
    end
end

%LLA figure
LLA_fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (17<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        LLA_fig_num = temp;
        flag_do_plots = 1;
    end
end






%% Find the grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize concatenation arrays
concatenate_scanLine_rings = [];
concatenate_LiDAR_XYZ_points = [];
concatenate_VehiclePose_XYZ_points = [];
concatenate_unit_ortho_vehicle_vectors_XYZ = [];
concatenate_LIDAR_intensity = [];



empty_XYZ = [nan nan nan];

% Loop through the scans we want to keep
for scanLineNumber = scanLineNumber_start:scanLineNumber_end
    LIDAR_scan = LiDAR_Scan_ENU_Entire_Loop{scanLineNumber};
    LIDAR_XYZ  = LIDAR_scan(:,1:3);
    LIDAR_intensity = LIDAR_scan(:,4);
    LIDAR_ringID = round(LIDAR_scan(:,5));
    
    indicies_to_keep = [];
    for ith_ring = 1:length(rings_to_analyze)
        this_ring_indicies = find(LIDAR_ringID==rings_to_analyze(ith_ring));
        indicies_to_keep = [indicies_to_keep; this_ring_indicies]; %#ok<AGROW>
    end
    Nindicies = length(indicies_to_keep);

    LIDAR_XYZ_to_add = [LIDAR_XYZ(indicies_to_keep,:); empty_XYZ];
    LIDAR_intensity_to_add = [LIDAR_intensity(indicies_to_keep,:); nan];
    scanLine_rings_to_add = [scanLineNumber*ones(Nindicies,1) LIDAR_ringID(indicies_to_keep,1); nan nan];
    Npoints_added = length(LIDAR_XYZ_to_add(:,1));


    concatenate_LiDAR_XYZ_points = [concatenate_LiDAR_XYZ_points; LIDAR_XYZ_to_add]; %#ok<AGROW>
    concatenate_LIDAR_intensity  = [concatenate_LIDAR_intensity; LIDAR_intensity_to_add]; %#ok<AGROW>
    concatenate_VehiclePose_XYZ_points = [concatenate_VehiclePose_XYZ_points; ones(Npoints_added,1)*VehiclePose(scanLineNumber,1:3)]; %#ok<AGROW>
    concatenate_unit_ortho_vehicle_vectors_XYZ = [...
        concatenate_unit_ortho_vehicle_vectors_XYZ; ...
        ones(Npoints_added,1)*[unit_ortho_vehicle_vectors_XY(scanLineNumber,:) 0]]; %#ok<AGROW>
    concatenate_scanLine_rings = [...
        concatenate_scanLine_rings; ...
        scanLine_rings_to_add]; %#ok<AGROW>
end



%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if flag_do_plots
    figure(fig_num);
    hold on;
    grid on;
    axis equal

    intensity_min = min(concatenate_LIDAR_intensity);
    intensity_max = max(concatenate_LIDAR_intensity);

    plot(concatenate_LiDAR_XYZ_points(:,1),concatenate_LiDAR_XYZ_points(:,2),'.','Color',[0 0 1],'MarkerSize',1);

    xlabel('East position [m]');
    ylabel('North position [m]');


    %% Plot the LIDAR in 3D ENU
figure(ENU_3D_fig_num);
clf;

hold on;
grid on;
axis equal

if 1==simple_view_toggle_1
    % Plot the vehicle's trajectory
    plot3(VehiclePose(:,1),VehiclePose(:,2),VehiclePose(:,3),'-','Color',[1 0 1],'MarkerSize',10,'LineWidth',3);

    % Plot start and end points of trajectory
    plot3(VehiclePose(1,1),VehiclePose(1,2),VehiclePose(1,3),'.','Color',[1 0 0],'MarkerSize',10,'LineWidth',3);
    plot3(VehiclePose(end,1),VehiclePose(end,2),VehiclePose(end,3),'o','Color',[0 1 0],'MarkerSize',10,'LineWidth',3);

else
    plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
    plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',3);

    % Show the orthogonal arrows showing vehicle motion directions. Green
    % is forward, bLue is Left
    quiver3(...
        VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
        unit_vehicle_change_in_pose_XY(scanLineNumber_start,1),unit_vehicle_change_in_pose_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 1 0]);
    quiver3(...
        VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
        unit_ortho_vehicle_vectors_XY(scanLineNumber_start,1),unit_ortho_vehicle_vectors_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
end

xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)

if 1==simple_view_toggle_1
    % Plot the LIDAR data simply as blue points
    plot3(concatenate_LiDAR_XYZ_points(:,1),concatenate_LiDAR_XYZ_points(:,2),concatenate_LiDAR_XYZ_points(:,3), '.','Color',[0 0 1],'MarkerSize',1);
else
    intensity_fraction = scaling_1*concatenate_LIDAR_intensity/(intensity_max - intensity_min);
    
    % Use user-defined colormap_string to map intensity to colors. For a
    % full example, see fcn_geometry_fillColorFromNumberOrName
    old_colormap = colormap;
    color_ordering = colormap('hot');
    colormap(old_colormap);
    N_colors = length(color_ordering(:,1));

    % Make sure the plot number is a fraction between 0 and 1
    plot_number = min(max(0,intensity_fraction),1);

    % Convert the plot number to a row
    color_row = floor((N_colors-1)*plot_number) + 1;


    % Plot the LIDAR data with intensity
    for ith_color = min(color_row):max(color_row)
        % Find the color
        color_vector = color_ordering(ith_color,:);

        % Find all the points that are in this color
        index_in_this_color = find(color_row==ith_color);
        plot3(...
            concatenate_LiDAR_XYZ_points(index_in_this_color,1),...
            concatenate_LiDAR_XYZ_points(index_in_this_color,2),...
            concatenate_LiDAR_XYZ_points(index_in_this_color,3), '.','Color',color_vector,'MarkerSize',5);
    end
end

set(gca,'CameraViewAngle',6)

%% Plot the LIDAR in XZ ENU
figure(ENU_XZ_fig_num);
clf;

hold on;
grid on;
axis equal

xlabel('East position [m]');
ylabel('Up position [m]');

% Plot the LIDAR data
plot(concatenate_LiDAR_XYZ_points(:,1),concatenate_LiDAR_XYZ_points(:,3), '.','Color',[0 0 1],'MarkerSize',5);

% plot the vehicle
plot(VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot(VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',1);

%% Plot the LIDAR in YZ ENU


gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class


figure(ENU_YZ_fig_num);
clf;

hold on;
grid on;
axis equal

xlabel('North position [m]');
ylabel('Up position [m]');

% Plot the LIDAR data
plot(concatenate_LiDAR_XYZ_points(:,2),concatenate_LiDAR_XYZ_points(:,3), '.','Color',[0 0 1],'MarkerSize',5);

% plot the vehicle
plot(VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot(VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',1);

% Plot the vehicle pose in LLA
figure(LLA_fig_num);
clf;

LLA_VehiclePose = gps_object.ENU2WGSLLA(VehiclePose(:,1:3));

% Do plot, and set it up so that zoom is correct, tick marks are correct,
% map is centered where we want, etc. To see options, do the geoplot, then
% zoom into the location we want, and then do:
%
% format long % This sets the format so we can see all the digits
% temp = gca  % This Gets Current Axis (GCA) and saves it as temp
%
% Next, click on "show all properties" and copy out the ZoomLevel and
% MapCenter values into the correct locations below. This way, every time
% we run this script, it automatically zooms and centers onto the correct
% location.

h_geoplot = geoplot(LLA_VehiclePose(:,1),LLA_VehiclePose(:,2),'-','Color',[0 0 1],'MarkerSize',10);
hold on;
h_parent =  get(h_geoplot,'Parent');
set(h_parent,'ZoomLevel',20.5,'MapCenter',[40.865718697633348 -77.830965127435817]);
geobasemap satellite
geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default


% NOTE: transition the geoplotting to use the PlotTestTrack library 
% Plot start and end points
geoplot(LLA_VehiclePose(1,1),LLA_VehiclePose(1,2),'.','Color',[0 1 0],'MarkerSize',10);
geoplot(LLA_VehiclePose(end,1),LLA_VehiclePose(end,2),'o','Color',[1 0 0],'MarkerSize',10);

% Plot the LIDAR in LLA
% Use the class to convert LLA to ENU
concatenate_LiDAR_LLA_points = gps_object.ENU2WGSLLA(concatenate_LiDAR_XYZ_points(:,1:3));


% Plot the LLA of LIDAR points
figure(LLA_fig_num);

if 1==simple_view_toggle_1
    % Plot the LIDAR data simply as magenta and black points
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'mo','MarkerSize',10);
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'k.','MarkerSize',10);
else
    intensity_fraction = scaling_2*concatenate_LIDAR_intensity/(intensity_max - intensity_min);
      
    % Use user-defined colormap_string to map intensity to colors. For a
    % full example, see fcn_geometry_fillColorFromNumberOrName
    old_colormap = colormap;
    % color_ordering = colormap('hot');
    color_ordering = flipud(colormap('sky'));
    colormap(old_colormap);
    N_colors = length(color_ordering(:,1));

    % Make sure the plot number is a fraction between 0 and 1
    plot_number = min(max(0,intensity_fraction),1);

    % Convert the plot number to a row
    color_row = floor((N_colors-1)*plot_number) + 1;


    % Plot the LIDAR data with intensity
    for ith_color = min(color_row):max(color_row)
        % Find the color
        color_vector = color_ordering(ith_color,:);

        % Find all the points that are in this color
        index_in_this_color = find(color_row==ith_color);

        geoplot(concatenate_LiDAR_LLA_points(index_in_this_color,1),concatenate_LiDAR_LLA_points(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',5);
    end
end

% Plot the vehicle pose on top of this
geoplot(LLA_VehiclePose(scanLineNumber_start:scanLineNumber_end,1),LLA_VehiclePose(scanLineNumber_start:scanLineNumber_end,2),'.','Color',[1 1 0],'MarkerSize',10);

%%

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

end