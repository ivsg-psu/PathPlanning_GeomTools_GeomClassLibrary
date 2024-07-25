% This code is used to transfer the output from
% script_test_geometry_updatedSurfacedAnalysis
% to the input of script_temp_Aneesh_example

fig_num = 1;
% Find the max and min of X and Y
max_x = max(true_boundary_points(:,1));
min_x = min(true_boundary_points(:,1));
max_y = max(true_boundary_points(:,2));
min_y = min(true_boundary_points(:,2));


% Find the length of x and y coordination, in this example, 
% 0.625 is the distance between points. 
x_range = (min_x:0.625:max_x)'; 
y_range = (min_y:0.625:max_y)';

% Find the indices of each point in term of X and Y
[~, indice_X] = ismember(true_boundary_points(:,1),x_range);
[~, indice_Y] = ismember(true_boundary_points(:,2),y_range);



z = zeros(length(x_range),length(y_range));

% Convert subscript indices to linear indices
linear_indices = sub2ind(size(z), indice_Y, indice_X);

% Assign the value 1 to the specified positions
z(linear_indices) = 1;

 
border_only_test_grid = z;
drive_path_rows_columns = gridCenters_driven_path;
true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);




% Dr B code 
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
end % Ends fcn_INTERNAL_findTrueBorders