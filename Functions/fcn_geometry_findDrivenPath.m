function [concatenate_LiDAR_XYZ_points_under_vehicle, concatenate_LiDAR_LLA_points_under_vehicle] = fcn_geometry_findDrivenPath(base_coordinates_LLA,... 
concatenate_LiDAR_XYZ_points, concatenate_VehiclePose_XYZ_points, concatenate_unit_ortho_vehicle_vectors_XYZ, vehicle_half_width, varargin)
%% fcn_geometry_findDrivenPath
%
% This function calculates the driven path and plots them on a given
% reference location. 
% 
% FORMAT: 
% [concatenate_LiDAR_XYZ_points_under_vehicle,
% concatenate_LiDAR_LLA_points_under_vehicle] =
% fcn_geometry_findDrivenPath(base_coordinates_LLA,
% concatenate_LiDAR_XYZ_points, concatenate_VehiclePose_XYZ_points,
% concatenate_unit_ortho_vehicle_vectors_XYZ, vehicle_half_width, varargin)
%
% INPUTS:
% base_coordinates_LLA
% concatenate_LiDAR_XYZ_points
% concatenate_VehiclePose_XYZ_points
% concatenate_unit_ortho_vehicle_vectors_XYZ
% vehicle_half_width
%
% (OPTIONAL INPUTS):
% flag_plot_in_3D
% fig_num_3D
% fig_num
%
%
% Revision history:
% 2024_07_17 - S Brennan
% -- wrote the code
% 2024_07_21 - Aleksandr 
% -- Functionalized this code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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
        narginchk(5,8);
 
    end
end

% Does user want to plot in 3D?
flag_plot_in_3D = 0;
if (6<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        flag_plot_in_3D = temp;
        if flag_plot_in_3D~= 0 && flag_plot_in_3D~= 1 
            error('The flag_plot_in_3D should be either 1 or 0')
        end
    end
end

fig_num_3D = [];
if (7<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        fig_num_3D = temp;
    end
end

% Does user want to specify fig_num?
LLA_fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (8<= nargin)
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

%% Find the driven path
% This must be done in ENU coordinates because LLA is not an orthogonal
% coordinate system. To find the surface, we find the distance of the lidar
% points in the orthogonal direction by taking the dot product of the LIDAR
% points, relative to the vehicle center with the unit projection vector
% pointed to the left of the vehicle.

% Calculate the vectors
vector_from_vehicle_pose_to_LIDAR_points = concatenate_LiDAR_XYZ_points - concatenate_VehiclePose_XYZ_points;

% Calculate the transverse distance
transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*concatenate_unit_ortho_vehicle_vectors_XYZ,2);

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
% lane_half_width = (3.6576/2) * 0.40;  
% vehicle_half_width = (3.6576/2) * 0.40;  

% Find the indices of the points under the vehicle.
indices_under_vehicle = find(abs(transverse_only_LIDAR_points)<vehicle_half_width);

% Find the points under the vehicle 
concatenate_LiDAR_XYZ_points_under_vehicle = concatenate_LiDAR_XYZ_points(indices_under_vehicle,:);

% Create a GPS object to convert ENU points to LLA
gps_object = GPS(base_coordinates_LLA(1),base_coordinates_LLA(2),base_coordinates_LLA(3)); % Load the GPS class

% Convert LiDAR ENU points under vehicle into LLA
concatenate_LiDAR_LLA_points_under_vehicle = gps_object.ENU2WGSLLA(concatenate_LiDAR_XYZ_points_under_vehicle(:,1:3));

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


    % Plot the LIDAR data underneath the vehicle in LLA
    figure(LLA_fig_num);
    geoplot(concatenate_LiDAR_LLA_points_under_vehicle(:,1),concatenate_LiDAR_LLA_points_under_vehicle(:,2),'g.','MarkerSize',2);
    hold on;
    grid on;

end % Ends check if plotting 

if flag_plot_in_3D

    % Plot the LIDAR data underneath the vehicle in XYZ
    figure(fig_num_3D)
    plot3(concatenate_LiDAR_XYZ_points_under_vehicle(:,1),concatenate_LiDAR_XYZ_points_under_vehicle(:,2),concatenate_LiDAR_XYZ_points_under_vehicle(:,3), '.','Color',[0 1 0],'MarkerSize',1);


    hold on;
    grid on;
    axis equal
    xlabel('East position [m]');
    ylabel('North position [m]');
    zlabel('Up position [m]');
    view(3)

end


if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end



end % Ends main function

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
