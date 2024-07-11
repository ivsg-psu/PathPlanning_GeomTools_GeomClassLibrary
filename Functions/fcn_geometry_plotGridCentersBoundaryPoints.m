function LLA_data = fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,varargin)
% This function takes in various boundary points and then generates a plot
% that shows where the boundary,driveable/non-driveable area, and unmapped area
% 
% FORMAT:  
% [LLA] =
% fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(fig_num))
%
% INPUTS:
% ENU_data: Data points 
% marker_size: size of the points that are going to be plotted
% RGB_triplet: Color of the markers
% (OPTIONAL INPUT)
% fig_num: figure number
%
% OUTPUTS:
% LLA coordinates of each point
%
% DEPENDENCIES:
% GPS CLASS
%
% EXAMPLES:
% See script: script_fcn_geometry_plotGridCentersBoundaryPoints
%
% Revision History
% 2024_07_10 - Aneesh Batchu 
% -- wrote the code originally
% 2024_07_10 - Aleksandr Goncharov
% -- Funtionlized this code
% 2024_07_11 - Aleksandr Goncharov
% -- Shortened and cleaned up the function to be more universal
%

flag_do_debug = 0; % Flag to plot the results for debugging
flag_max_speed = 0; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; 
else
    debug_fig_num = []; 
end
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

if flag_max_speed==0
    % Are there the right number of inputs?
    narginchk(3,4);
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (4<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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

% % True boundary points
% Trace_coordinates = [true_boundary_points,zeros(length(true_boundary_points),1)]; 

 % Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_data = gps_object.ENU2WGSLLA(ENU_data);

% % Unmapped grid centers with zero point density in LLA
% LLA_gridCenters_zero_point_density = gps_object.ENU2WGSLLA(gridCenters_zero_point_density);
% 
% % Unmapped grid centers with low point density in LLA
% LLA_gridCenters_low_point_density = gps_object.ENU2WGSLLA(gridCenters_low_point_density);
% 
% % Mapped grid centers 
% LLA_gridCenters_mapped_grids = gps_object.ENU2WGSLLA(gridCenters_mapped_grids(:,1:3));
% 
% % Drivable grid centers
% drivable_grid_centers_ENU = gridCenters_mapped_grids((gridCenters_mapped_grids(:,4) == 1),1:3); 
% LLA_gridCenters_drivable_grids = gps_object.ENU2WGSLLA(drivable_grid_centers_ENU);
% 
% % Non-drivable grid centers
% non_drivable_grid_centers_ENU = gridCenters_mapped_grids((gridCenters_mapped_grids(:,4) == 0),1:3); 
% LLA_gridCenters_non_drivable_grids = gps_object.ENU2WGSLLA(non_drivable_grid_centers_ENU);
% 
% LLA=[LLA_data_computed_boundary_pts,LLA_gridCenters_zero_point_density,LLA_gridCenters_low_point_density,LLA_gridCenters_drivable_grids,LLA_gridCenters_non_drivable_grids];

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

% Plot the LLA boundary points
figure(fig_num);
clf;
hold on
xlabel('Latitude')
ylabel('Longitude')
hold off 

% Plot
geoplot(LLA_data(:,1),LLA_data(:,2),'.','MarkerSize',marker_size,'Color',RGB_triplet);

% % Plot the unmapped grid centers with low point density
% geoplot(LLA_gridCenters_low_point_density(:,1),LLA_gridCenters_low_point_density(:,2),'.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
% 
% % Plot the mapped grid centers 
% geoplot(LLA_gridCenters_mapped_grids(:,1),LLA_gridCenters_mapped_grids(:,2),'.','MarkerSize',30,'Color',[0.3 0.3 0.3]);
% 
% % Plot the mapped grid centers 
% geoplot(LLA_gridCenters_drivable_grids(:,1),LLA_gridCenters_drivable_grids(:,2),'g.','MarkerSize',15);
% 
% % Plot the mapped grid centers 
% geoplot(LLA_gridCenters_non_drivable_grids(:,1),LLA_gridCenters_non_drivable_grids(:,2),'r.','MarkerSize',15);
% 
% % Plot the hand-labeled boundary points
% geoplot(hand_labeled_boundary_points_LLA(:,1),hand_labeled_boundary_points_LLA(:,2),'y.','MarkerSize',40);
% % % hold on
% % geoplot(boundary_points_LLA(:,1),boundary_points_LLA(:,2),'k.','MarkerSize',10);
% 
% % Plot the computed boundary points
% % geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'y.','MarkerSize',40);
% geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'c.','MarkerSize',30);
% geoplot(LLA_data_computed_boundary_pts(:,1),LLA_data_computed_boundary_pts(:,2),'b.','MarkerSize',15);

title('Boundary Points in LLA ')
%hold off
geobasemap satellite
geotickformat -dd

end
end
