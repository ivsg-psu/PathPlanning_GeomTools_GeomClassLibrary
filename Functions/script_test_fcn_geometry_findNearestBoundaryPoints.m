% script_test_fcn_geometry_findNearestBoundaryPoints
% Exercises the function: fcn_geometry_findNearestBoundaryPoints
% Revision history:
% 2024_7_25
% Jiabao Zhao wrote the code


%% Test 1, simple example, three points
gridCenters_driven_path = [3 1;2 2;1 3];
gridCenters_non_drivable_grids = [1 1;3 3];
grid_size = 1;
true_boundary_points = [1 1; 3 3];
fig_num = 1;
[true_borders,x_borders,y_borders] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
    gridCenters_non_drivable_grids, gridCenters_driven_path, grid_size,fig_num);
assert(isequal(length(true_borders(:,1)),4));
assert(isequal(length(x_borders),4));
assert(isequal(length(y_borders),4));


%% Test 2  Real Data 

fig_num = 1224; 
% after running script_test_geometry_updatedSurfaceAnalysis
[true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
    gridCenters_non_drivable_grids, gridCenters_driven_path, fig_num);

%%

[~, C] = ismember(true_boundary_points(:,1),x_range);



%% DATA

gridCenters_driven_path = [-112.3080   55.5933
 -112.3080   56.3933
 -111.5080   56.3933
 -112.3080   57.1933
 -111.5080   57.1933
 -111.5080   57.9933
 -110.7080   57.9933
 -111.5080   58.7933
 -110.7080   58.7933
 -110.7080   59.5933
 -109.9080   59.5933
 -110.7080   60.3933
 -109.9080   60.3933
 -109.9080   61.1933
 -109.1080   61.1933
 -109.1080   61.9933
 -108.3080   61.9933
 -109.1080   62.7933
 -108.3080   62.7933
];


gridCenters_non_drivable_grids = [-109.1080   53.9933         0
 -109.1080   54.7933         0
 -114.7080   55.5933         0
 -109.1080   55.5933         0
 -108.3080   55.5933         0
 -113.9080   56.3933         0
 -108.3080   56.3933         0
 -113.9080   57.1933         0
 -113.1080   57.1933         0
 -108.3080   57.1933         0
 -107.5080   57.1933         0
 -113.1080   57.9933         0
 -107.5080   57.9933         0
 -113.1080   58.7933         0
 -112.3080   58.7933         0
 -107.5080   58.7933         0
 -106.7080   58.7933         0
 -112.3080   59.5933         0
 -106.7080   59.5933         0
 -112.3080   60.3933         0
 -111.5080   60.3933         0
 -105.9080   60.3933         0
 -111.5080   61.1933         0
 -105.9080   61.1933         0
 -111.5080   61.9933         0
 -110.7080   61.9933         0
 -110.7080   62.7933         0
 -109.9080   63.5933         0
];
true_boundary_points = [-113.1080   56.7933
 -112.3080   58.3933
 -111.5080   59.9933
 -110.7080   61.5933
 -109.9080   63.1933
 -109.1080   55.9933
 -108.3080   57.5933
 -107.5080   59.1933
 -106.7080   59.9933
 -109.5080   54.7933
 -109.5080   55.5933
 -108.7080   56.3933
 -108.7080   57.1933
 -107.9080   57.9933
 -107.9080   58.7933
 -107.1080   59.5933
 -106.3080   60.3933
 -106.3080   61.1933
 -113.5080   56.3933
 -112.7080   57.1933
 -112.7080   57.9933
 -111.9080   58.7933
 -111.9080   59.5933
 -111.1080   60.3933
 -111.1080   61.1933
 -110.3080   61.9933
 -110.3080   62.7933];

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

% Find the indices of each point in term of X and Y range
[~,indice_X] = ismember(true_boundary_points(:,1),x_range);
[~,indice_Y] = ismember(true_boundary_points(:,2),y_range);
z = zeros(length(x_range),length(y_range));

% Return the indice of corresponding points with 1
z(sub2ind(size(z), indice_Y, indice_X)) = 1;
border_only_test_grid = z;



% Find the max and min of X and Y coordination
max_x = max(gridCenters_non_drivable_grids(:,1));
min_x = min(gridCenters_non_drivable_grids(:,1));
max_y = max(gridCenters_non_drivable_grids(:,2));
min_y = min(gridCenters_non_drivable_grids(:,2));

% Find the range of x and y coordination and break them into grids
x_range = (min_x:(grid_size):max_x)'; 
y_range = (min_y:(grid_size):max_y)';

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


% Find the indices of each point in term of X and Y range
[~, indice_X] = ismember(new_gridCenters_driven_path(:,1),x_range);
[~, indice_Y] = ismember(new_gridCenters_driven_path(:,2),y_range);
drive_path_rows_columns = [indice_Y,indice_X];

true_borders_indices = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);
true_borders_indices_unique = unique(true_borders_indices,'rows','stable');
true_borders_x = x_range(true_borders_indices_unique(:,1));
true_borders_y = y_range(true_borders_indices_unique(:,2));
true_borders = [true_borders_x true_borders_y];



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
    pause_duration = 0.01;
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
        pause(pause_duration);
    end

    % Loop up row-wise
    for ith_row = current_row:Nrows
        if flag_do_debug
            set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
            pause(pause_duration);
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
            pause(pause_duration);
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
            pause(pause_duration);
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
            pause(pause_duration);
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
end % Ends fcn_INTERNA
