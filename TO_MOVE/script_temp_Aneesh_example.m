
%% Test 1: a simple corner
%%%%
% Set up an empty grid
grid_width = 20;
grid_height = 8;
empty_test_grid = zeros(grid_height,grid_width);
Y_values = (1:grid_height)';
X_values = (1:grid_width)';
[column_subscripts, row_subscripts] = meshgrid(X_values,Y_values);

% Flip the rows up/down so that the matrix looks the same as XY
% row_subscripts = flipud(row_subscripts);
row_indicies = reshape(row_subscripts,[],1);
column_indicies = reshape(column_subscripts,[],1);


%%%%
% Set an artifial boundary
% Set all data to boundaries above y = 2*x+3
indicies_as_borders = (row_indicies>=(2*column_indicies));
border_only_test_grid = empty_test_grid;
border_only_test_grid(indicies_as_borders) = 1;


%%%%
% Set a driving area
% Set the x values
drive_path_column = (2:grid_width)';

% Force the path to follow y = 0.5*z, and limit it to range of grid height
% and 0:
drive_path_rows = round(0.5*drive_path_column);
bad_data = drive_path_rows>grid_height;
drive_path_rows(bad_data) = [];
drive_path_column(bad_data) = [];
drive_path_rows_columns = [drive_path_rows drive_path_column];

%%%%
% Find the true borders
true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);


%% Test 2: an isolated square with the corner
%%%%
% Set up an empty grid
grid_width = 20;
grid_height = 8;
empty_test_grid = zeros(grid_height,grid_width);
Y_values = (1:grid_height)';
X_values = (1:grid_width)';
[column_subscripts, row_subscripts] = meshgrid(X_values,Y_values);

% Flip the rows up/down so that the matrix looks the same as XY
% row_subscripts = flipud(row_subscripts);
row_indicies = reshape(row_subscripts,[],1);
column_indicies = reshape(column_subscripts,[],1);


%%%%
% Set an artifial boundary
% Set all data in a specific square to 1
indicies_as_borders = (row_indicies>=(2*column_indicies));
border_only_test_grid = empty_test_grid;
border_only_test_grid(indicies_as_borders) = 1;
border_only_test_grid(1:3,15:18) = 1;


%%%%
% Set a driving area
% Set the x values
drive_path_column = (2:grid_width)';

% Force the path to follow y = 0.5*z, and limit it to range of grid height
% and 0:
drive_path_rows = round(0.5*drive_path_column);
bad_data = drive_path_rows>grid_height;
drive_path_rows(bad_data) = [];
drive_path_column(bad_data) = [];
drive_path_rows_columns = [drive_path_rows drive_path_column];

%%%%
% Find the true borders
true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);

%% Test 3: a blocking square
%%%%
% Set up an empty grid
grid_width = 20;
grid_height = 8;
empty_test_grid = zeros(grid_height,grid_width);
Y_values = (1:grid_height)';
X_values = (1:grid_width)';
[column_subscripts, row_subscripts] = meshgrid(X_values,Y_values);

% Flip the rows up/down so that the matrix looks the same as XY
% row_subscripts = flipud(row_subscripts);
row_indicies = reshape(row_subscripts,[],1);
column_indicies = reshape(column_subscripts,[],1);


%%%%
% Set an artifial boundary
% Set all data in a specific square to 1
indicies_as_borders = (row_indicies>=(2*column_indicies));
border_only_test_grid = empty_test_grid;
border_only_test_grid(indicies_as_borders) = 1;
border_only_test_grid(1:3,15:18) = 1;
border_only_test_grid(6:7,6:8) = 1;

%%%%
% Set a driving area
% Set the x values
drive_path_column = (2:grid_width)';

% Force the path to follow y = 0.5*z, and limit it to range of grid height
% and 0:
drive_path_rows = round(0.5*drive_path_column);
bad_data = drive_path_rows>grid_height;
drive_path_rows(bad_data) = [];
drive_path_column(bad_data) = [];
drive_path_rows_columns = [drive_path_rows drive_path_column];

%%%%
% Find the true borders
true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);

%% Test 1: a simple corner
%%%%
% Set up an empty grid
grid_width = 20;
grid_height = 8;

start_num = 2; 

empty_test_grid = zeros(grid_height,grid_width);
Y_values = (start_num:(grid_height+(start_num-1)))';
X_values = (start_num:(grid_width+(start_num-1)))';
[column_subscripts, row_subscripts] = meshgrid(X_values,Y_values);

% Flip the rows up/down so that the matrix looks the same as XY
% row_subscripts = flipud(row_subscripts);
row_indicies = reshape(row_subscripts,[],1);
column_indicies = reshape(column_subscripts,[],1);


%%%%
% Set an artifial boundary
% Set all data to boundaries above y = 2*x+3
indicies_as_borders = (row_indicies>=(2*column_indicies));
border_only_test_grid = empty_test_grid;
border_only_test_grid(indicies_as_borders) = 1;


%%%%
% Set a driving area
% Set the x values

drive_path_start_column = 2; 

drive_path_column = (drive_path_start_column:grid_width)';

% Force the path to follow y = 0.5*z, and limit it to range of grid height
% and 0:
drive_path_rows = round(0.5*drive_path_column);
bad_data = drive_path_rows>grid_height;
drive_path_rows(bad_data) = [];
drive_path_column(bad_data) = [];
drive_path_rows_columns = [drive_path_rows drive_path_column];

%%%%
% Find the true borders
true_borders = fcn_INTERNAL_findTrueBorders(border_only_test_grid,drive_path_rows_columns, fig_num);


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

%% fcn_INTERNAL_findTrueBorders
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