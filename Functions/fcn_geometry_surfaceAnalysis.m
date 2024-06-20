function [drivable_grids, non_drivable_grids, unmapped_grids] = fcn_geometry_surfaceAnalysis(inputPoints, gridSize, gridBoundaries, pointDensity, varargin)


% EXAMPLES:
%
% See the script: script_test_fcn_geometry_surfaceAnalysis
% for a full test suite.

% This function was written on 2024_06_17 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu 

% Revision history:
% 2024_01_15 - Aneesh Batchu
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(4,5);
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

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
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

flag_do_debug = 1; 

% Divides the data into grids
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

% Calculate the number of points in each grid
total_N_points_in_each_grid = fcn_geometry_findRepeatedIndices(gridIndices,length(grid_AABBs(:,1)), -1);

% Save all the indices to a cell array. Each cell contains the indices of
% the inputPoints that belong to the grid. 
gridIndices_cell_array = fcn_geometry_createGridPointIndicesCellArray(gridIndices,length(grid_AABBs(:,1)),-1);

% Unmapped/not_fitted grids. These grids do not contain enough number of
% points to fit a plane and classify it as drivable or non-drivable
original_unmapped_grid_numbers = find(total_N_points_in_each_grid(:,1) < pointDensity); 

% Output
unmapped_grids = original_unmapped_grid_numbers; 

% Mapped grids. Later these grids are classified into drivable and non
% drivable
original_mapped_grid_numbers = find(total_N_points_in_each_grid(:,1) >= pointDensity); 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_grids = gridIndices_cell_array(original_mapped_grid_numbers); 

% Indices of points in mapped grids
indices_original_mapped_grids = vertcat(original_mapped_grids{:}); 

% Input points in the mapped grids
points_in_original_mapped_grids = inputPoints(indices_original_mapped_grids,:); 


if flag_do_debug

    plot_3D = 0; 

    fig_num = 100011; 
    figure(fig_num); clf;

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

    % The indices of the mapped grids are extracted and concatenated
    original_unmapped_grids = gridIndices_cell_array(original_unmapped_grid_numbers);

    % Indices of points in mapped grids
    indices_original_unmapped_grids = vertcat(original_unmapped_grids{:});

    % Input points in the mapped grids
    points_in_original_unmapped_grids = inputPoints(indices_original_unmapped_grids,:);
 
    if plot_3D
        view(3)
        % Plot the unmapped points red
        plot3(points_in_original_unmapped_grids(:,1),points_in_original_unmapped_grids(:,2),points_in_original_unmapped_grids(:,3),'.','MarkerSize',20,'Color',[0.6350 0.0780 0.1840]);
        % Plot the mapped points green
        plot3(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),points_in_original_mapped_grids(:,3),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);
    else
        % Plot the unmapped points red
        plot(points_in_original_unmapped_grids(:,1),points_in_original_unmapped_grids(:,2),'.','MarkerSize',20,'Color',[0.6350 0.0780 0.1840]);
        % Plot the mapped points green
        plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);
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
end

% Total number of mapped grids
total_mapped_grids = length(original_mapped_grids); 

% Parameters of the plane fits of each mapped grid
parameters_of_fitted_plane = zeros(total_mapped_grids,3); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 

% A loop that iterates based on the length of mapped grids and fits the
% plane to each mapped grid
for ith_mapped_grid = 1:total_mapped_grids
    [parameters, standard_deviation_in_z, z_fit, unit_normal_vectors(ith_mapped_grid,:), base_point, standard_deviation_in_plane_orthogonals(ith_mapped_grid,:)] = fcn_geometry_fitPlaneLinearRegression(inputPoints(original_mapped_grids{ith_mapped_grid},:),-1);
        parameters_of_fitted_plane(ith_mapped_grid,:) = parameters';
end

% Plane analysis: setting thresholds

% STEP 1
% Comparing normal vector with verticle direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2); 

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product); 

% Find the grids (with a fitted plane) that are within the vertical
% threshold (change to a different name: vertical threshold)
mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < pi/9;

% STEP 2: Standard deviation of the orthogonal (perpendicular) distances of
% the points to the plane (after fit) 
% Find the grids that are within standard deviation limit
% This is not enough (delta Y) is also important
mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < 0.1; 

% Grids that satisy both the conditions (STEP 1 and STEP 2). The grids that
% are within the vertical and std threshold
mapped_grids_within_vertical_and_std_thresholds = (mapped_grids_within_vertical_threshold == 1) & (mapped_grids_within_std_threshold == 1);

% Find the drivable grids (original)
drivable_grids = original_mapped_grid_numbers(mapped_grids_within_vertical_and_std_thresholds); 

% Find the non-drivable grids (original)
non_drivable_grids = original_mapped_grid_numbers(mapped_grids_within_vertical_and_std_thresholds == 0);

% Final drivable grid numbers of the mapped grids
drivable_grid_numbers_in_mapped_grids = find(ismember(original_mapped_grid_numbers, drivable_grids));

% Final non drivable grid numbers of the mapped grids
non_drivable_grid_numbers_in_mapped_grids = find(ismember(original_mapped_grid_numbers, non_drivable_grids));

if flag_do_debug
    fig_num = 100012;
    figure(fig_num); clf;

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

    view(3)

%     % Plot all the points
% plot(inputPoints(:,1),inputPoints(:,2),'.','MarkerSize',20,'Color',[0 0 0]);
% 

    % Plot all the points
    plot3(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),points_in_original_mapped_grids(:,3),'.','MarkerSize',20,'Color',[0 0 0]);

    for ith_domain = 1:length(drivable_grids)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
        current_color = [0.4660 0.6740 0.1880];

        % Plot current AABB
        current_AABB = grid_AABBs(drivable_grids(ith_domain),:);

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
        rows_in_domain = gridIndices==drivable_grids(ith_domain);
        points_in_domain = inputPoints(rows_in_domain,:);
        plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);

    % Plot the unit normal vector
    quiver3(gridCenters(drivable_grids(ith_domain),1),gridCenters(drivable_grids(ith_domain),2),gridCenters(drivable_grids(ith_domain),3), unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),1),unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),2),unit_normal_vectors(drivable_grid_numbers_in_mapped_grids(ith_domain),3),0,'g','Linewidth',3);
    
    end

    for ith_domain = 1:length(non_drivable_grids)
        % Get current color
        % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
        current_color = [0.6350 0.0780 0.1840];

        % Plot current AABB
        current_AABB = grid_AABBs(non_drivable_grids(ith_domain),:);

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
        rows_in_domain = gridIndices==non_drivable_grids(ith_domain);
        points_in_domain = inputPoints(rows_in_domain,:);
        plot3(points_in_domain(:,1),points_in_domain(:,2),points_in_domain(:,3),'.','MarkerSize',20,'Color',current_color);

    % Plot the unit vector
    quiver3(gridCenters(non_drivable_grids(ith_domain),1),gridCenters(non_drivable_grids(ith_domain),2),gridCenters(non_drivable_grids(ith_domain),3), unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),1),unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),2),unit_normal_vectors(non_drivable_grid_numbers_in_mapped_grids(ith_domain),3),0,'g','Linewidth',3);
    end

    % % Plot the unit vector
    % quiver3(base_point(1,1),base_point(1,2),base_point(1,3), unit_vector(1,1),unit_vector(1,2),unit_vector(1,3),0,'g','Linewidth',3);

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
% if flag_do_plots
% 
% 
%     temp_h = figure(fig_num);
%     flag_rescale_axis = 0;
%     if isempty(get(temp_h,'Children'))
%         flag_rescale_axis = 1;
%     end
% 
%     hold on;
%     grid on;
%     axis equal
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     zlabel('Z [m]')
% 
% 
%     % Make axis slightly larger?
%     if flag_rescale_axis
%         temp = axis;
%         %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
%         axis_range_x = temp(2)-temp(1);
%         axis_range_y = temp(4)-temp(3);
%         percent_larger = 0.3;
%         axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
%     end
% 
% 
% end % Ends check if plotting

% if flag_do_debug
%     fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
% end

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