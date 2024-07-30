function [ ...
        current_grid_numbers_of_driven_path, ...
        total_points_in_each_grid_in_the_driven_path, ...
        total_points_in_each_grid_with_points_greater_than_zero, ...
        gridCenters_driven_path ...
    ] = fcn_geometry_plotGridCentersAndBoundaryPoints (...
        gridCenters_greater_than_zero_point_density, ...
        grids_greater_than_zero_points, ...
        boundary_points_driven_path, ...
        fig_num_gridCenters_and_boundary_points_greater_than_zero ...
    )
%% Debugging and Input checks

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


%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "inpolygon" is used to find the grids within the boundary points 
[in,~] = inpolygon(gridCenters_greater_than_zero_point_density(:,1),gridCenters_greater_than_zero_point_density(:,2), ...
                   boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = grids_greater_than_zero_points(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = find(in); 

% Total points in each grid in the driven pa
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(grids_greater_than_zero_points); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2)];


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
