function [true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points, ...
     gridCenters_driven_path, varargin)
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
%      true_borders_x: the x coordiantion of true_borders
%
%      true_borders_y: the y coordiantion of true_borders
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
flag_max_speed = 0;
if (nargin== 3&& isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(1,3);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 2<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

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
max_x = max(max(true_boundary_points(:,1)),max(gridCenters_driven_path(:,1)));
min_x = min(min(true_boundary_points(:,1)),min(gridCenters_driven_path(:,1)));
max_y = max(max(true_boundary_points(:,2)),max(gridCenters_driven_path(:,2)));
min_y = min(min(true_boundary_points(:,2)),min(gridCenters_driven_path(:,2)));

% Find the range of x and y coordination and break them into grids
% Round the number to four decimal places so the function could read every
% decimals.
x_range = (min_x:grid_size:max_x)';
y_range = (min_y:grid_size:max_y)';
x_rounded_true_boundary_points = round((true_boundary_points(:,1)), 4);
rounded_x_range = round((x_range), 4);
y_rounded_true_boundary_points = round((true_boundary_points(:,2)), 4);
rounded_y_range = round((y_range), 4);

% Find the indices of each point in term of X and Y range of the boundary
% points
[~,indice_X] = ismember(x_rounded_true_boundary_points,rounded_x_range);
[~,indice_Y] = ismember(y_rounded_true_boundary_points,rounded_y_range);

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

% Find the indices of each point in term of X and Y range of driving path
[~, indice_X1] = ismember(round_new_gridCenters_driven_path(:,1),rounded_x_range);
[~, indice_Y1] = ismember(round_new_gridCenters_driven_path(:,2),rounded_y_range);

% remove any rows that contain zeros for indices
drive_path_rows_columns = [indice_Y1,indice_X1];
rows_to_keep = all(drive_path_rows_columns~=0,2);
% one of the input
drive_path_rows_columns = drive_path_rows_columns(rows_to_keep, :);

% Add a boundary point at the starting indices of the driving path if
% needed
if indice_X ~= indice_X1(1)
    added_point = [indice_X1(1) length(y_range)];
    indice_X = [indice_X ; added_point(1,1)];
    indice_Y = [indice_Y ; added_point(1,2)];
end

% make a zeros matrix that is large enough to fit all the data
z = zeros(max(max(length(indice_Y)),max(length(y_range))),max(max(length(indice_X)),max(length(x_range))));
% Return the indice of corresponding points with 1
z(sub2ind(size(z), indice_Y, indice_X)) = 1;

% one of the input
border_only_test_grid = z;
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
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    Nrows = length(border_only_test_grid(:,1));
    Ncols = length(border_only_test_grid(1,:));


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

    % Again, need to be very careful here - the "path" is in XY coordinates and
    % is plotted as such. So columns are "X" and rows are "Y".

    % Loop through all the path points
    true_borders_indices = [];

    if flag_do_debug
        h_test_drive_path = plot(nan,nan,'o','Color',[1 0 0],'MarkerSize',20);
        h_test_expansion  = plot(nan,nan,'.','Color',[0 0 1],'MarkerSize',20);
        %pause_duration = 0.01;
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
                true_borders_indices = [true_borders_indices; ith_row current_column]; %#ok<AGROW>
                break
            end
        end

        if flag_do_debug
            plot(true_borders_indices(:,2), true_borders_indices(:,1),'g.','MarkerSize',10);
        end

        % Loop down row-wise
        for ith_row = current_row:(-1):1
            if flag_do_debug
                set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
                %pause(pause_duration);
            end
            if 1==border_only_test_grid(ith_row,current_column)
                true_borders_indices = [true_borders_indices; ith_row current_column]; %#ok<AGROW>
                break
            end
        end

        if flag_do_debug
            plot(true_borders_indices(:,2), true_borders_indices(:,1),'g.','MarkerSize',10);
        end

        % Loop up column-wise
        for jth_column = current_column:Ncols
            if flag_do_debug
                set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
                %  pause(pause_duration);
            end
            if 1==border_only_test_grid(current_row, jth_column)
                true_borders_indices = [true_borders_indices; current_row jth_column]; %#ok<AGROW>
                break
            end
        end

        if flag_do_debug
            plot(true_borders_indices(:,2), true_borders_indices(:,1),'g.','MarkerSize',10);
        end

        % Loop down column-wise
        for jth_column = current_column:(-1):1
            if flag_do_debug
                set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
                %pause(pause_duration);
            end
            if 1==border_only_test_grid(current_row, jth_column)
                true_borders_indices = [true_borders_indices; current_row jth_column]; %#ok<AGROW>
                break
            end
        end

        if flag_do_debug
            plot(true_borders_indices(:,2), true_borders_indices(:,1),'g.','MarkerSize',10);
        end

    end
    if ~isempty(fig_num)
        % Plot the results
        plot(true_borders_indices(:,2),true_borders_indices(:,1),'o','Color',[0 1 0],'MarkerSize',20)
    end
    % find the true boundary points
    true_borders_indices_unique = unique(true_borders_indices,'rows','stable');
    true_borders_x = y_range(true_borders_indices_unique(:,1));
    true_borders_y = x_range(true_borders_indices_unique(:,2));
    true_borders = [true_borders_x true_borders_y];
    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
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
