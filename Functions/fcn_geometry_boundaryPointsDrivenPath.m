function [boundary_points_driven_path,shift_distance]...
    = fcn_geometry_boundaryPointsDrivenPath( ...
        VehiclePose, ...
        varargin)
% This function is written to find the boundary points of driven path
%
% FORMAT: 
%
%       [z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference]...
%       = fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map...
%       ,(file_name)(fig1),(fig2),(fig3));
%
% INPUTS:
%
%       cell_start: starting cell
%       cell_end: ending cell
%       color_map: color map to plot the points in
%       (OPTIONAL INPUTS)
%
%       file_name: if an external file is needed for LiDAR_ENU_cell
%       fig_1: Label for figure 1
%       fig_2: Label for figure 2
%       fig_3: Label for figure 3
%
% OUTPUTS:
%
%       z, mean_z, deviation_from_mean
%       angles, mean_angle,angle_difference
%
% DEPENDENCIES:
%
%       The script must either have existing LiDAR_ENU_cell or import a file
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_plotLidarDeviation
%
%
% Revision History
% Aneesh Batchu - 2024_07_18
% -- wrote the code originally
% Steven Young 2024_07_24
% -- functionalized the code
% 




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

if flag_max_speed == 0
    % Are there the right number of inputs?
    narginchk(3,4);
end

%file name if needed
if 4<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        file_name=temp;%#ok<NASGU>
    end
end

%labeling figures

flag_do_plots = 0;
fig_1=[];
fig_2=[];
fig_3=[];

%figure 1
if flag_max_speed==0
if 5<= nargin
    temp = varargin{2};
    if ~isempty(temp)
        fig_1=temp;
        flag_do_plots=1;
    end
end

%figure 2
if 6<= nargin
    temp = varargin{3};
    if ~isempty(temp)
        fig_2=temp;
        flag_do_plots=1;
    end
end
%figure 3
if 7 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_3=temp;
        flag_do_plots=1;
    end
end
end
%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

% Transverse distance of the left and right boundary points of the driven
% path from vehicle center 
transverse_distance_of_boundary_points = [lane_half_width*unit_ortho_vehicle_vectors_XY, zeros(length(unit_ortho_vehicle_vectors_XY),1)];

shift = 2.5; 
% Shift
shift_distance = [unit_vehicle_change_in_pose_XY*shift, zeros(length(unit_vehicle_change_in_pose_XY),1)]; 

% % Left boundary points of the driven path
% left_boundary_points = VehiclePose(:,1:3) + transverse_distance_of_boundary_points - shift_distance; 
% 
% % Left boundary points of the driven path
% right_boundary_points = VehiclePose(:,1:3) - transverse_distance_of_boundary_points - shift_distance; 

% Left boundary points of the driven path
left_boundary_points = VehiclePose(:,1:3) + transverse_distance_of_boundary_points; 

% Left boundary points of the driven path
right_boundary_points = VehiclePose(:,1:3) - transverse_distance_of_boundary_points; 



boundaryLineNumber_start = 1400 - 8; 
boundaryLineNumber_end = 1410 - 6; 

% lengths_boundary_points = [lane_half_width*unit_ortho_vehicle_vectors_XY(scanLineNumber_start:scanLineNumber_end,:), zeros((scanLineNumber_end-scanLineNumber_start+1),1)];
% 
% left_boundary_points = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:3) + lengths_boundary_points;
% right_boundary_points = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:3) - lengths_boundary_points;

fig_num = 1098;
figure(fig_num);clf;

hold on;
grid on;
axis equal

plot3(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot3(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',3);

% Show the orthogonal arrows showing vehicle motion directions. Green
% is forward, bLue is Left
quiver3(...
    VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
    unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,1),unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 1 0]);
quiver3(...
    VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
    unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,1),unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);

plot3(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 1 0],'MarkerSize',30,'LineWidth',3);
plot3(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 1],'MarkerSize',30,'LineWidth',3);


xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)


%%

ENU_3D_fig_num = 3;
figure(ENU_3D_fig_num);

plot3(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 1 0],'MarkerSize',30,'LineWidth',3);
plot3(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 1],'MarkerSize',30,'LineWidth',3);

boundary_points_driven_path = [right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3);
    flipud(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3));
    right_boundary_points(boundaryLineNumber_start,1:3)];

boundary_points_driven_path_LLA = gps_object.ENU2WGSLLA(boundary_points_driven_path);

LLA_fig_num = 2;
figure(LLA_fig_num);
% Plot the LIDAR data underneath the vehicle in LLA
figure(LLA_fig_num);
geoplot(boundary_points_driven_path_LLA(:,1),boundary_points_driven_path_LLA(:,2),'b.','MarkerSize',30);


% %%  INPOLYGON
% boundary_points = [right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3); flipud(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3)); right_boundary_points(boundaryLineNumber_start,1:3)];
% 
% points = VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3);
% 
% % xv = [0.5;0.2;1.0;0;0.8;0.5];
% % yv = [1.0;0.1;0.7;0.7;0.1;1];
% % 
% % 
% % xq = [0.1;0.5;0.9;0.2;0.4;0.5;0.5;0.9;0.6;0.8;0.7;0.2];
% % yq = [0.4;0.6;0.9;0.7;0.3;0.8;0.2;0.4;0.4;0.6;0.2;0.6];
% 
% 
% [in,on] = inpolygon(points(:,1),points(:,2),boundary_points(:,1),boundary_points(:,2));
% 
% figure(122)
% 
% plot(boundary_points(:,1),boundary_points(:,2)) % polygon
% 
% hold on
% plot(points(in,1),points(in,2),'r+') % points strictly inside
% % plot(xq(on),yq(on),'k*') % points on edge
% plot(points(~in,1),points(~in,2),'bo') % points outside
% hold off
%% Any debugging?
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
    if ~isempty(fig_1)
    figure(fig_1)
    % Create the scatter plot
    scatter3(x, y, z, 20, deviation_from_mean, 'filled');
    colormap(color_map);
    colorbar;
    view(2);
    title('Deviation of z-value from the mean');
    xlabel('x');
    ylabel('y');
    end
    
    if ~isempty(fig_2)
    figure(fig_2)
    % Create the scatter plot
    scatter3(x, y, z, 20, angle_difference, 'filled');
    colormap(color_map);
    colorbar;
    view(2);
    title('Angle difference from the mean angle');
    xlabel('x');
    ylabel('y');
    end
    
    if ~isempty(fig_3)
    figure(fig_3)
    colors = z;
    scatter3(x,y,z,20,colors);
    colormap(color_map); % You can use 'viridis' if you have it, or other colormaps
    colorbar; % Add a color bar to show the color mapping
    view(2)
    end
    
end
end