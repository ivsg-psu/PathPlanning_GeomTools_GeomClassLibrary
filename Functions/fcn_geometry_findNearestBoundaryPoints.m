function [true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points, ...
     gridCenters_driven_path, fig_num)
%% fcn_geometry_findNearestBoundaryPoints
% Find the nearest boundary points of a road 
% 
% FORMAT:
%
%      [true_borders] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,
%      gridCenters_non_drivable_grids, gridCenters_driven_path, (fig _num))
%
% INPUTS:     
%       
%      true_boundary_points: The boundary points between the non-drivable 
%      path and the drivable path.
%
%      gridCenters_non_drivable_grids: Grid points that are non-drivable
%      for vehicles.
%
%      gridCenters_driven_path: Grid points that are drivable for vehicles.
%
%      gird_size: size of the grid     
%
%      (OPTIONAL INPUT)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%       
%      true_borders: the nearest boundary points that are from the center of 
%      driven path.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_findNearestBoundaryPoints.m
%       test suite.
%
% This function was written on 2024_07_25 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.


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

% Find the grid size
Low_to_high = sort(true_boundary_points(:,1));
unique_low_to_high = unique(Low_to_high);
grid_size = abs(unique_low_to_high(2)-unique_low_to_high(1));

% Find the max and min of X and Y coordination 
max_x = max(true_boundary_points(:,1));
min_x = min(true_boundary_points(:,1));
max_y = max(true_boundary_points(:,2));
min_y = min(true_boundary_points(:,2));


% Find the range of x and y coordination and break them into grids
x_range = (min_x:grid_size:max_x)'; 
y_range = (min_y:grid_size:max_y)';
x_rounded_true_boundary_points = round((true_boundary_points(:,1)), 4);
rounded_x_range = round((x_range), 4);
y_rounded_true_boundary_points = round((true_boundary_points(:,2)), 4);
rounded_y_range = round((y_range), 4);

% Find the indices of each point in term of X and Y range
[~,indice_X] = ismember(x_rounded_true_boundary_points,rounded_x_range);
[~,indice_Y] = ismember(y_rounded_true_boundary_points,rounded_y_range);
z = zeros(length(indice_Y),length(indice_X));
% Return the indice of corresponding points with 1
z(sub2ind(size(z), indice_Y, indice_X)) = 1;
border_only_test_grid = z;

% Offset distance of a driven grid from a 
offset_distance = grid_size;

% Compute new points directly using matrix operations
% Initialize matrices for new points
top_points_driven_path = gridCenters_driven_path + [0, offset_distance];
bottom_points_driven_path = gridCenters_driven_path - [0, offset_distance];
left_points_driven_path = gridCenters_driven_path - [offset_distance, 0];
right_points_driven_path = gridCenters_driven_path + [offset_distance, 0];
new_gridCenters_driven_path = [gridCenters_driven_path;top_points_driven_path;
    bottom_points_driven_path; right_points_driven_path;left_points_driven_path];

round_new_gridCenters_driven_path = round(new_gridCenters_driven_path, 4);

% Find the indices of each point in term of X and Y range
[~, indice_X] = ismember(round_new_gridCenters_driven_path(:,1),rounded_x_range);
[~, indice_Y] = ismember(round_new_gridCenters_driven_path(:,2),rounded_y_range);
drive_path_rows_columns = [indice_Y,indice_X];

true_borders_indices = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);

% find the true boundary points
true_borders_indices_unique = unique(true_borders_indices,'rows','stable');
true_borders_x = y_range(true_borders_indices_unique(:,1));
true_borders_y = x_range(true_borders_indices_unique(:,2));
true_borders = [true_borders_x true_borders_y];
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

%% %% fcn_INTERNAL_findTrueBorders 
function true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num)
flag_do_debug = 1;

Nrows = length(border_only_test_grid(:,1));
Ncols = length(border_only_test_grid(1,:));

if flag_do_debug && ~isempty(fig_num)
    %%%
    % Plot the inputs
    figure(38383);
    clf;
    hold on;

    % Plot all the empty grids as points. Note: rows are "height" and columns
    % are "width". Need to be VERY careful keeping track of this difference.
    matrix_indicies = reshape(border_only_test_grid,[],1);
    matrix_indicies_to_plot = matrix_indicies==0;
   
    row_indicies = repmat((1:Nrows)',Ncols,1);
    column_indicies = reshape(repmat((1:Ncols),Nrows,1),[],1);

    plot(column_indicies(matrix_indicies_to_plot),row_indicies(matrix_indicies_to_plot), '.','Color',[0.8 0.8 0.8]);
    xlim([0 Ncols+1])
    ylim([0 Nrows+1])
    xlabel('Columns (X)');
    ylabel('Rows (Y)');

    % Plot the non-empty values as black dots
    matrix_indicies = reshape(border_only_test_grid,[],1);
    matrix_indicies_to_plot = matrix_indicies==1;
    plot(column_indicies(matrix_indicies_to_plot),row_indicies(matrix_indicies_to_plot),'.','Color',[0 0 0],'MarkerSize',20)

    % Plot the results. Again, be careful: columns are X, rows are Y
    plot(drive_path_rows_columns(:,2),drive_path_rows_columns(:,1),'.','Color',[1 0 0],'MarkerSize',20)

end

% Again, need to be very careful here - the "path" is in XY coordinates and
% is plotted as such. So columns are "X" and rows are "Y".

% Loop through all the path points
true_borders = [];

if flag_do_debug
    h_test_drive_path = plot(nan,nan,'o','Color',[1 0 0],'MarkerSize',20);
    h_test_expansion  = plot(nan,nan,'.','Color',[0 0 1],'MarkerSize',20);
    % pause_duration = 0.01;
end

for ith_drive_path_point = 1:length(drive_path_rows_columns)

    % For each path point, find row and column of the current point. All
    % the loops will start at this point and search "outward" for
    % true_borders.
    current_point_to_test = drive_path_rows_columns(ith_drive_path_point,:);
    current_row = current_point_to_test(1,1);
    current_column = current_point_to_test(1,2);
    if flag_do_debug
        set(h_test_drive_path,'Xdata',current_column,'Ydata',current_row);
        % pause(pause_duration);
    end

    % Loop up row-wise
    for ith_row = current_row:Nrows
        if flag_do_debug
            set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
            % pause(pause_duration);
        end
        if 1==border_only_test_grid(ith_row,current_column)
            true_borders = [true_borders; ith_row current_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug
        plot(true_borders(:,2), true_borders(:,1),'g.','MarkerSize',10);
    end

    % Loop down row-wise
    for ith_row = current_row:(-1):1
        if flag_do_debug
            set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
            % pause(pause_duration);
        end
        if 1==border_only_test_grid(ith_row,current_column)
            true_borders = [true_borders; ith_row current_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug
        plot(true_borders(:,2), true_borders(:,1),'g.','MarkerSize',10);
    end

    % Loop up column-wise
    for jth_column = current_column:Ncols
        if flag_do_debug
            set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
            % pause(pause_duration);
        end
        if 1==border_only_test_grid(current_row, jth_column)
            true_borders = [true_borders; current_row jth_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug
        plot(true_borders(:,2), true_borders(:,1),'g.','MarkerSize',10);
    end

    % Loop down column-wise
    for jth_column = current_column:(-1):1
        if flag_do_debug
            set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
            % pause(pause_duration);
        end
        if 1==border_only_test_grid(current_row, jth_column)
            true_borders = [true_borders; current_row jth_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug
        plot(true_borders(:,2), true_borders(:,1),'g.','MarkerSize',10);
    end

end
if ~isempty(fig_num)
    % Plot the results
    plot(true_borders(:,2),true_borders(:,1),'o','Color',[0 1 0],'MarkerSize',20)
end
end % Ends fcn_INTERNAL_findTrueBorders