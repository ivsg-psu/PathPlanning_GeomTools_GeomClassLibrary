function [gridIndicies, grid_AABBs, gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries,varargin)
% Introduction:
% This function aims to seprate given points in X/XY/XYZ format into a
% user defined grid size and grid boundaries and output a cell array
% containing indices of points in the domain,the domain in AABB format and
% the grid centers
% Example: a set of XYZ points, grid space of two meters, the function will
% seprate all the points in [0,2m], [2m,4m], etc. If the points are in 2D
% there should be an option in the function to create the grid in x-axis or
% in y-axis
%
% FORMAT:
%
%       [gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num))
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%
%     none
%
% OUTPUTS:
%
%     gridIndicies : a matrix of size [N x 1] constaining the domain number
%     of the point. If a point is not in any domain, it is NaN.
%
%     grid_AABBs : matrix in the format of AABB (i.e axis aligned bounding
%     box format) that gives the domains corresponding to the indices of
%     points present in the domain. The format of the AABB is 
%
%             [x_low x_high y_low y_high z_low z_high] 
% 
%     for each row, with one row per domain.
%
%     gridCenters : matrix in 2D or 3D that gives the respective
%     gridDomain centers
%
%     lefover_points: points outside grid boundaries
%
%
% DEPENDENCIES:
%
%      fcn_geometry_fillColorFromNumberOrName
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_separatePointsIntoGrids.m for a full
%       test suite.
%
% This function was written on 2024_01_22 by V. Wagh, updated by S. Brennan
% Questions or comments? vbw5054@psu.edu or sbrennan@psu.edu


% Revision history:
% 2024_01_22 by V. Wagh (vbw5054@psu.edu)
% -- start writing code
% 2024_01_24 by V. Wagh
% -- changed inputs to include gridBoundaries
% -- changed outputs
% 2024_01_26 - S. Brennan
% -- added 2D plotting
% -- cleaned up comments a bit
% -- added fast mode
% -- added optional figure number
% -- added plotting for all cases
% -- added support for 1D
% -- fixed name on script and function
% 2024_04_14 - S. Brennan
% -- added fcn_geometry_fillColorFromNumberOrName

flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
    debug_fig_num = [];  %#ok<NASGU> % no plotting yet
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
        narginchk(3,4); % no optional inputs


        % Check the inputPoints input

        % Test if is numeric input
        if ~isnumeric(inputPoints)
            error('The inputPoints must be numeric')
        end
        size_of_points = size(inputPoints);
        if size_of_points(2)<1 || size_of_points(2)>3
            error('Gridding only allowed for 1D, 2D, or 3D points.')
        end
        
        % Check the gridSize input
        fcn_DebugTools_checkInputsToFunctions(gridSize, 'positive_1column_of_numbers',1);

         
        % Check the size_of_vector input

          % Test if is numeric input
        if ~isnumeric(gridBoundaries)
            error('The gridBoundaries must be numeric')
        end
        size_of_vector = size(gridBoundaries);
        if ~isequal(size_of_vector,[1 2*size_of_points(2)])
            error('The gridBoundaries must have dimension of [1 x 2*N] where N is the dimension, in format of: [low_x high_x low_y high_y low_z high_z]')
        end
    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if (0==flag_max_speed)
    if (4 == nargin)
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
            figure(fig_num);
            flag_do_plots = 1;
        end
    elseif flag_do_debug
        fig = figure;
        fig_num = fig.Number;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many dimensions do we have?
N_dims = round(length(gridBoundaries(1,:))/2);


% Calculate the grid boundaries from the gridBoundaries input
min_max_form = reshape(gridBoundaries,2,N_dims)';

% Calculate the number of grids in each dimension
nGrids = ceil((min_max_form(:,2) - min_max_form(:,1))/gridSize);

% total number of domains i.e total number of grids
totalDomains = prod(nGrids,1);

% Find the end points of each dimension
start_end_grid_boundaries = [min_max_form(:,1) min_max_form(:,1)+gridSize*(nGrids-1)];

% Find the number of times each matrix needs to be repeated. This is
% related to the cumulative product of the counts.
cumprod_up   = cumprod(nGrids);
cumprod_down = cumprod(flipud(nGrids));
row_repeats = [1; cumprod_up(1:end-1)];
col_repeats = [flipud(cumprod_down(1:end-1)); 1];

% Permute the elements along each grid direction. The way this operates is
% it examines each row of grid boundaries, in each dimension, to calculate
% how many times it has to repeat within a row form, and within a column
% form, to fill out the permutation matrix.
grid_AABBs  = zeros(totalDomains,N_dims*2);
gridCenters  = zeros(totalDomains,N_dims*2);
for ith_dimension = 1:N_dims

    % What is the vector to repeat in this dimension?
    vector_to_freeze = (start_end_grid_boundaries(ith_dimension,1):gridSize:start_end_grid_boundaries(ith_dimension,2));

    % Repeat it by the requisite number of rows
    vector_repeated_by_rows = repmat(vector_to_freeze,row_repeats(ith_dimension),1);

    % Reshape it into a column
    n_elements = numel(vector_repeated_by_rows);
    vector_reshaped_as_column = reshape(vector_repeated_by_rows,n_elements,1);

    % Repeat the column into a super-column
    vector_repeated_by_columns = repmat(vector_reshaped_as_column,col_repeats(ith_dimension),1);

    % Save the results
    grid_AABBs(:,ith_dimension*2-1) = vector_repeated_by_columns;
    grid_AABBs(:,ith_dimension*2)   = vector_repeated_by_columns+gridSize;
    gridCenters(:,ith_dimension)     = vector_repeated_by_columns + gridSize/2;

end


% calculating indices:
gridIndicies = nan(length(inputPoints(:,1)),1);
good_point_rows = find(prod((inputPoints>=min_max_form(:,1)').*(inputPoints<=min_max_form(:,2)'),2));

min_point = min_max_form(:,1)';
gridSubscripts = floor((inputPoints(good_point_rows,:)-min_point)./gridSize) + 1;

if N_dims==1
    gridIndicies(good_point_rows) = gridSubscripts(:,1);
elseif N_dims==2
    gridIndicies(good_point_rows) = sub2ind(nGrids',gridSubscripts(:,1),gridSubscripts(:,2));
elseif N_dims==3
    gridIndicies(good_point_rows) = sub2ind(nGrids',gridSubscripts(:,1),gridSubscripts(:,2),gridSubscripts(:,3));
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


    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    % Plot the grid with the points
    if N_dims==1 % 1D points

        % Plot all the points
        plot(inputPoints(:,1),inputPoints(:,1)*0,'.','MarkerSize',20,'Color',[0 0 0]);

        % Plot the input points by domain with different colors for each
        % domain
        for ith_domain = 1:totalDomains
            % Get current color
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

            % Plot current AABB
            current_AABB = grid_AABBs(ith_domain,:);

            % Nudge the current AABB inward
            current_AABB = current_AABB + gridSize/100*[1 -1];

            % Calculate the gridlines
            gridlines = [...
                current_AABB(1,1); ...
                current_AABB(1,2); ...
                nan;
                ];

            % Plot the result
            plot(gridlines(:,1),0*gridlines(:,1),'-','Color',current_color,'LineWidth',3);

            % Get all points in this domain and plot them
            rows_in_domain = gridIndicies==ith_domain;
            points_in_domain = inputPoints(rows_in_domain,:);
            plot(points_in_domain(:,1),0*points_in_domain(:,1),'.','MarkerSize',20,'Color',current_color);
        end

    elseif N_dims==2 % 2D points

        % Plot all the points
        plot(inputPoints(:,1),inputPoints(:,2),'.','MarkerSize',20,'Color',[0 0 0]);

        % Plot the input points by domain with different colors for each
        % domain
        for ith_domain = 1:totalDomains
            % Get current color
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

            % Plot current AABB
            current_AABB = grid_AABBs(ith_domain,:);

            % Nudge the current AABB inward
            current_AABB = current_AABB + gridSize/100*[1 -1 1 -1];

            % Calculate the gridlines
            gridlines = [...
                current_AABB(1,1) current_AABB(1,3); ...
                current_AABB(1,1) current_AABB(1,4); ...
                nan nan;
                current_AABB(1,2) current_AABB(1,3); ...
                current_AABB(1,2) current_AABB(1,4); ...
                nan nan;
                current_AABB(1,1) current_AABB(1,3); ...
                current_AABB(1,2) current_AABB(1,3); ...
                nan nan;
                current_AABB(1,1) current_AABB(1,4); ...
                current_AABB(1,2) current_AABB(1,4); ...
                ];

            % Plot the result
            plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

            % Get all points in this domain and plot them
            rows_in_domain = gridIndicies==ith_domain;
            points_in_domain = inputPoints(rows_in_domain,:);
            plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',20,'Color',current_color);
        end



    else % 3D vectors
        zlabel('Z [m]')

        view(3);

        % Plot all the points
        plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

        % Plot the input points by domain with different colors for each
        % domain
        for ith_domain = 1:totalDomains
            % Get current color
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

            % Plot current AABB
            current_AABB = grid_AABBs(ith_domain,:);

            % Nudge the current AABB inward
            current_AABB = current_AABB + gridSize/100*[1 -1 1 -1 1 -1];

            % Calculate the gridlines
            gridlines = [...
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
                nan nan nan;
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,5); ...
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,5); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,5); ...
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,5); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,3) current_AABB(1,6); ...
                current_AABB(1,2) current_AABB(1,3) current_AABB(1,6); ...
                nan nan nan;
                current_AABB(1,1) current_AABB(1,4) current_AABB(1,6); ...
                current_AABB(1,2) current_AABB(1,4) current_AABB(1,6); ...
                nan nan nan];

            % Plot the result
            plot3(gridlines(:,1),gridlines(:,2),gridlines(:,3),'-','Color',current_color,'LineWidth',3);


            % Get all points in this domain and plot them
            rows_in_domain = gridIndicies==ith_domain;
            points_in_domain = inputPoints(rows_in_domain,:);
            plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);
        end

    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    
end % Ends check if plotting

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


%% old code
%     % Initialize cell array to store grid data
%     gridData = cell(numGridsX, numGridsY);
%
%     % Iterate over each grid
%     for i = 1:numGridsX
%         for j = 1:numGridsY
%             % Define the boundaries of the current grid
%             gridXMin = xMin + (i - 1) * gridSize;
%             gridXMax = xMin + i * gridSize;
%             gridYMin = yMin + (j - 1) * gridSize;
%             gridYMax = yMin + j * gridSize;
%
%             % Find indices of points within the current grid
%             pointsInGrid = find(inputPoints(:, 1) >= gridXMin & inputPoints(:, 1) < gridXMax & ...
%                                 inputPoints(:, 2) >= gridYMin & inputPoints(:, 2) < gridYMax);
%
%             % Store grid data in the cell array
%             gridData{i, j}.indices = pointsInGrid;
%             gridData{i, j}.boundaries = [gridXMin, gridXMax, gridYMin, gridYMax];
%             gridData{i, j}.center = [mean([gridXMin, gridXMax]), mean([gridYMin, gridYMax])];
%         end
%     end
% end
